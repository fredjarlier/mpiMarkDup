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
     log.h
     
   Authors:
    Frederic Jarlier,   Institut Curie
    Firmain Martin,     Paris Descartes University
*/

#ifndef MD_LOG_H
#define MD_LOG_H

/**
 * @file log.h
 */

#include <mpi.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#define LOG_USE_COLOR

#define LOG_BUF_SIZ 4096

#define MPI_ASSERT(fx) do { \
  assert(fx == MPI_SUCCESS); \
} while(0)

/// Log levels.
enum mdLogLevel {
    MD_LOG_OFF = -1,  ///< All logging disabled.
    MD_LOG_ERROR,     ///< Logging of errors only.
    MD_LOG_WARNING,   ///< Logging of errors and warnings.
    MD_LOG_INFO,      ///< Logging of errors, warnings, and normal but significant events.
    MD_LOG_DEBUG,     ///< Logging of all except the most detailed debug events.
    MD_LOG_TRACE      ///< All logging enabled.
};

#ifdef LOG_USE_COLOR
static const char *level_colors[] = {
    "\x1b[35m", "\x1b[31m", "\x1b[33m", "\x1b[32m", "\x1b[36m", "\x1b[94m"
};
#endif

/// Sets the selected log level.
void md_set_log_level(enum mdLogLevel level);
void md_set_log_comm(MPI_Comm comm) ;

/// Gets the selected log level.
enum mdLogLevel md_get_log_level();
MPI_Comm md_get_log_comm() ;

/// Selected log level.
/**
 * One of the MD_LOG_* values. The default is MD_LOG_WARNING.
 * @note Avoid direct use of this variable. Use md_set_log_level and md_get_log_level instead.
 */
extern int md_verbose;

void md_log_rank(int rank, enum mdLogLevel severity, const char *file,  const char *context, const int line,  const char *format, ...) ;
void md_log_all(enum mdLogLevel severity, const char *file,  const char *context, const int line, const char *format, ...) ;

#define md_log_error(format, ...)   \
 if ( MD_LOG_ERROR  <= md_get_log_level()) {   \
    md_log_rank(0, MD_LOG_ERROR, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_warning(format, ...) \
 if ( MD_LOG_WARNING  <= md_get_log_level()) { \
    md_log_rank(0, MD_LOG_WARNING, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_info(format, ...)    \
 if ( MD_LOG_INFO  <= md_get_log_level()) {   \
    md_log_rank(0, MD_LOG_INFO, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_debug(format, ...)   \
 if ( MD_LOG_DEBUG <= md_get_log_level())  {  \
    md_log_rank(0, MD_LOG_DEBUG, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }
#define md_log_trace(format, ...)   \
 if ( MD_LOG_TRACE  <= md_get_log_level()) {  \
    md_log_rank(0, MD_LOG_TRACE, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__); \
 }

#define md_log_rank_error(rank, format, ...)   md_log_rank(rank, MD_LOG_ERROR, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_warning(rank, format, ...) md_log_rank(rank, MD_LOG_WARNING, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_info(rank, format, ...)    md_log_rank(rank, MD_LOG_INFO, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_debug(rank, format, ...)   md_log_rank(rank, MD_LOG_DEBUG, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_rank_trace(rank, format, ...)   md_log_rank(rank, MD_LOG_TRACE, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)

#define md_log_all_error(format, ...)          md_log_all(MD_LOG_ERROR, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_warning(format, ...)        md_log_all(MD_LOG_WARNING, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_info(format, ...)           md_log_all(MD_LOG_INFO, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_debug(format, ...)          md_log_all(MD_LOG_DEBUG, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)
#define md_log_all_trace(format, ...)          md_log_all(MD_LOG_TRACE, __FILE__, __func__, __LINE__, format, ##__VA_ARGS__)

#endif /* ifndef MD_LOG_H */
