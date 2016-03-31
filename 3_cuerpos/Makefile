##### Crea un pdf con las figuras
run : fig9.dat fig12.dat fig15.dat fig16.dat fig17.dat
	python figuras.py
##### Borra archivos creados
clean :
	rm *.x *.dat *.png
##### Datos Figura 9
9.x : fig9.c
	gcc $< -lm -o $@
fig9.dat : 9.x
	./$< > $@ 
##### Datos Figura 12
12.x : fig12.c
	gcc $< -lm -o $@
fig12.dat : 12.x
	./$< > $@

##### Datos Figura 15
15.x : fig15.c
	gcc $< -lm -o $@
fig15.dat : 15.x
	./$< > $@

##### Datos Figura 16
16.x : fig16.c
	gcc $< -lm -o $@
fig16.dat : 16.x
	./$< > $@

##### Datos Figura 17
17.x : fig17.c
	gcc $< -lm -o $@
fig17.dat : 17.x
	./$< > $@
