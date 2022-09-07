default: dif_rad2 dit_rad2 fht

dif_rad2: dif_rad2.c
	gcc dif_rad2.c -o dif_rad2 -lm -fopenmp

dit_rad2: dit_rad2.c
	gcc -std=c99 -Wall -g -lm -o dit_rad2 dit_rad2.c -fopenmp

fht: fht.c
	gcc -std=c99 -Wall -W -g -lm -o fht fht.c -fopenmp

clean:
	rm -f fht dit_rad2 dif_rad2
