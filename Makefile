cops:cops.c
	cc cops.c binhash.c interact.c state.c params.c -o cops -lm

clean:
	-rm cops