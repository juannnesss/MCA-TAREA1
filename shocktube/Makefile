
#Creacion solucion analitica y graficas

*.png : *.txt analitica.py
	python analitica.py

#Resultados de solucion numerica
*.txt: a.out
	./a.out

#Compilacion archivos en C

a.out: shocktube.c steps.c riemann.c
	gcc shocktube.c -lm
run: a.out *.txt *.png 
	
clean:
	rm -f *.txt *.png a.out
