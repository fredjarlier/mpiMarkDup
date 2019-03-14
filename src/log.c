/*
   mpiSORT
   Copyright (C) 2016-2019 Institut Curie / Institut Pasteur
   mpiSORT is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   mpiSORT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser Public License
   along with mpiSORT.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
   Module:
     log.c
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/


#include "log.h"

/**
 * @file log.c
 * @note inspired from https://github.com/pysam-developers/pysam/blob/master/htslib/htslib/hts_log.h
 * @todo use constructor to create communicator restricted logger.
 */

/**
 * logger severity level, shouldn't be use directly.
 * Use md_set_log_level/md_get_log_level instead.
 */


struct mdLog {
    int verbose;
    MPI_Comm comm;
};

struct mdLog md = {.verbose = MD_LOG_DEBUG, .comm = MPI_COMM_WORLD};

/**
 * @date 2018 Mar 24
 * @brief set the logger's level
 * @param[in] level logger severity level
 * @note Could be used to debug a snippet by setting md_verbose to MD_LOG_OFF,
 *       and call
 *       md_set_log_level(MD_LOG_TRACE);
 *           .... (snippet and log) ....
 *       md_set_log_level(MD_LOG_OFF);
 */

void md_set_log_level(enum mdLogLevel level) {
    md.verbose = level;
}

/**
 * @date 2018 Mar 24
 * @brief get the logger's level
 * @return logger severity level
 */

enum mdLogLevel md_get_log_level() {
    return md.verbose;
}

/**
 * @date 2018 Apr 2
 * @brief set the logger's communicator
 * @param[in] comm logger communicator
 */

void md_set_log_comm(MPI_Comm comm) {
    md.comm = comm;
}


/**
 * @date 2018 Mar 24
 * @brief get the logger's communicator
 * @return logger's communicator
 */

MPI_Comm md_get_log_comm() {
    return md.comm;
}

/**
 * @date 2018 Mar 24
 * @brief Given logger level, return its corresponding string
 * @param[in] level logger severity level
 * @return level corresponding string
 */

const char *get_level_tag(enum mdLogLevel level) {
    switch (level) {

        case MD_LOG_ERROR:
            return " ERROR ";

        case MD_LOG_WARNING:
            return " WARN  ";

        case MD_LOG_INFO:
            return " INFO  ";

        case MD_LOG_DEBUG:
            return " DEBUG ";

        case MD_LOG_TRACE:
            return " TRACE ";

        default:
            break;
    }

    return "*";
}

/**
 * @date                2018 Mar 24
 * @brief Logs event through all process.
 * @param level         Severity of the event:
 *                          - MD_LOG_OFF do nothing.
 *                          - MD_LOG_ERROR means that something went wrong so that a task could not be completed.
 *                          - MD_LOG_WARNING means that something unexpected happened, but that execution can continue, perhaps in a degraded mode.
 *                          - MD_LOG_INFO means that something normal but significant happened.
 *                          - MD_LOG_DEBUG means that something normal and insignificant happened.
 *                          - MD_LOG_TRACE means that something happened that might be of interest when troubleshooting.
 * @param file          use to pass __FILE__
 * @param function      use to pass __func__
 * @param line          use to pass __LINE__
 * @param format        Format string with placeholders, like printf.
 * @note                MPI's part inspired from https://stackoverflow.com/a/5310506
 * @note                If this function is not called by all process in the same snippet, we have DEADLOCK.
 * @todo                use buf_len to manage too long message, troncate it.
 * @todo                test if it's possible to communicate with a md_log_rank function
 */

void md_log_all(enum mdLogLevel level, const char *file,  const char *function, const int line,  const char *format, ...) {

    if (level <= md.verbose) {
        int myrank, num_procs;
        MPI_Comm_rank(md.comm, &myrank);
        MPI_Comm_size(md.comm, &num_procs);
        MPI_Status  status;

        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        va_list argptr;
        va_start(argptr, format);

        char buffer1[LOG_BUF_SIZ];
        char buffer2[LOG_BUF_SIZ];

        /* message's header */
#ifdef LOG_USE_COLOR

        if (level >= MD_LOG_TRACE) {
            sprintf(buffer1, "rank %4d :: [%02d:%02d:%02d][%s%s\x1b[0m] at %s:%s:%d :: ",  myrank, tm.tm_hour, tm.tm_min, tm.tm_sec, level_colors[level],
                    get_level_tag(level), file, function, line);

        } else {
            sprintf(buffer1, "rank %4d :: [%02d:%02d:%02d][%s%s\x1b[0m] ", myrank, tm.tm_hour, tm.tm_min, tm.tm_sec, level_colors[level], get_level_tag(level));
        }

#else

        if (level >= MD_LOG_TRACE) {
            sprintf(buffer1, "rank %4d :: [%02d:%02d:%02d][%s] at %s:%s:%d :: ", myrank, tm.tm_hour, tm.tm_min, tm.tm_sec, get_level_tag(level), file, function,
                    line);

        } else {
            sprintf(buffer1, "rank %4d :: [%02d:%02d:%02d][%s] ", myrank, tm.tm_hour, tm.tm_min, tm.tm_sec, get_level_tag(level));
        }

#endif
        vsprintf(buffer2, format, argptr);
        strcat(buffer1, buffer2);

        /* MPI's part, print in proc's order */

        if (myrank == 0) {
            int buf_len[num_procs];
            //C99, must *alloc if we need a portable version
            char msg[num_procs][LOG_BUF_SIZ];

            /* first print proc 0 msg */
            fprintf(stderr, "%s", buffer1);

            /* then print other proc's msg */
            for (int i = 1; i < num_procs; i++) {
                MPI_Recv(&buf_len[i], 1, MPI_INT, i, MPI_ANY_TAG, md.comm, &status);
                MPI_Recv(msg[i], LOG_BUF_SIZ, MPI_CHAR, i, MPI_ANY_TAG, md.comm, &status);
                fprintf(stderr, "%s", msg[i]);
            }

        } else {
            int buf_len = strlen(buffer1);
            MPI_Send(&buf_len, 1, MPI_INT, 0, myrank, md.comm);
            MPI_Send(buffer1, LOG_BUF_SIZ, MPI_CHAR, 0, myrank, md.comm);

        }
    va_end(argptr);
    }


}

/**
 * @date                2018 Mar 24
 * @brief Logs event through one process
 * @param rank          procs rank
 * @param level         Severity of the event:
 *                          - MD_LOG_OFF do nothing.
 *                          - MD_LOG_ERROR means that something went wrong so that a task could not be completed.
 *                          - MD_LOG_WARNING means that something unexpected happened, but that execution can continue, perhaps in a degraded mode.
 *                          - MD_LOG_INFO means that something normal but significant happened.
 *                          - MD_LOG_DEBUG means that something normal and insignificant happened.
 *                          - MD_LOG_TRACE means that something happened that might be of interest when troubleshooting.
 * @param file          use to pass __FILE__
 * @param function      use to pass __func__
 * @param line          use to pass __LINE__
 * @param format        Format string with placeholders, like printf.
 * @todo                exit when it is called for error level
 */

void md_log_rank(int rank, enum mdLogLevel level, const char *file,  const char *context, const int line,  const char *format, ...) {

    if (level <= md.verbose) {
        int myrank;
        MPI_Comm_rank(md.comm, &myrank);
        time_t t = time(NULL);
        struct tm tm = *localtime(&t);
        va_list argptr;
        va_start(argptr, format);
        char buffer[LOG_BUF_SIZ] = {0};

        if (myrank == rank) {
#ifdef LOG_USE_COLOR

            if (level >= MD_LOG_TRACE) {
                sprintf(buffer, "rank %4d :: [%02d:%02d:%02d][%s%s\x1b[0m] at %s:%s:%d :: ", rank, tm.tm_hour, tm.tm_min, tm.tm_sec, level_colors[level], get_level_tag(level), file, context, line);

            } else {
                sprintf(buffer, "rank %4d :: [%02d:%02d:%02d][%s%s\x1b[0m] ", rank, tm.tm_hour, tm.tm_min, tm.tm_sec, level_colors[level], get_level_tag(level));
            }

#else

            if (level >= MD_LOG_TRACE) {
                sprintf(buffer, "rank %4d :: [%02d:%02d:%02d][%s] at %s:%s:%d :: ", rank, tm.tm_hour, tm.tm_min, tm.tm_sec, get_level_tag(level), file, context, line);

            } else {
                sprintf(buffer, "rank %4d :: [%02d:%02d:%02d][%s] ", rank, tm.tm_hour, tm.tm_min, tm.tm_sec, get_level_tag(level));
            }

#endif
            strcat(buffer, format);
            vfprintf(stderr, buffer, argptr);
        }
    va_end(argptr);
    }

}
