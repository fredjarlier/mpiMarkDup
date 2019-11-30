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
     qksort.c

   Authors:
    Frederic Jarlier, 	Institut Curie
	Nicolas Joly, 		Institut Pasteur
	Nicolas Fedy,		Institut Curie
	Leonor Sirotti,	 	Institut Curie
	Thomas Magalhaes,	Institut Curie
	Paul Paganiban,		Institut Curie
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "qksort.h"

//size_t *base_arr2;

int compare_size_t(size_t *a, size_t *b, size_t *base_arr2){

	 size_t aa = *(size_t *)a, bb = *(size_t *)b;

	 if (base_arr2[aa] > base_arr2[bb])
		return 1;
	else if (base_arr2[aa] < base_arr2[bb])
		return -1;
	else
		return 0;

}

int compare_size_t_V2(const void *a, const void *b){

	if (*(const size_t *)a > *(const size_t *)b)
		return 1;
	else if (*(const size_t *)a < *(const size_t *)b)
		return -1;
	else
		return 0;

}

int partition(void *data, size_t *base_arr2, size_t i, size_t k, 
	int (*compare)(size_t *key1, size_t *key2, size_t *base_arr2)){

	size_t *a = data;
	size_t *pval, *temp;
	size_t r[3];
	/*
	 * allocate value for the partition value and swapping
	 */
	if ((pval = (size_t *)malloc(sizeof(size_t))) == NULL)
		return -1;
	if ((temp = (size_t *)malloc(sizeof(size_t))) == NULL){
		free(pval);
		return -1;
	}

	/*
	 * Use the median-of-three method to find partition value
	 */

	int part = (k - i + 1);
	if(part <= 0)
		part = 1;

	r[0] = (rand() % part +i);
	r[1] = (rand() % part +i);
	r[2] = (rand() % part +i);


	/*
	 * TODO: replace the qsort with issort
	 *
	 * issort(r, 3, sizeof(size_t), compare_size_t_V2);
	 */
	qsort(r, 3, sizeof(size_t),compare_size_t_V2);
	memcpy(pval, &a[r[1]], sizeof(size_t));

	/*
	 * Create 2 partitions around the partition value
	 */
	i--;
	k++;
	while(1) {

		/*
		 * move left until an element is found in the wrong partition
		 */

		do {
			k--;
		} while (compare(&a[k], pval, base_arr2) > 0);

		/*
		 * move right until an element is found in the wrong partition
		 */

		do {
			i++;
		} while (compare(&a[i], pval, base_arr2) < 0);

		if (i >= k){
			/*
			 * break when left and right counter cross
			 */
			break;
		}

		else{
			// swap element under the left and right counters
			memcpy(temp, &a[i], sizeof(size_t));
			memcpy(&a[i], &a[k], sizeof(size_t));
			memcpy(&a[k], temp, sizeof(size_t));
		}
	}

	/*
	 * free the storage allocated for partitioning
	 */
	free(pval);
	free(temp);

	/*
	 * return position dividing the two partition
	 */
	return k;

}

int qksort(void *data, size_t *base_arr2, size_t size, size_t i, size_t k, 
	         int (*compare)(size_t *key1, size_t *key2, size_t *base_arr2)){

	int j;

	/*
	 * stop recursion when it is not possible to partition further
	 * when calling qksort:
	 * i = 0
	 * k = size-1
	 */

	while (i < k){

		/*
		 * find where to partition the elements
		 */


		if ((j = partition(data, base_arr2, i, k, compare)) < 0){
			return -1;
		}
		/*
		 * recursively sort the left partition
		 */
		if (qksort(data, base_arr2, size, i, j, compare) < 0)
			return -1;
		/*
		 * iterate and sort the right partition
		 */
		i = j + 1;
	}
	return 0;
}
