compile:
	cd shocktube && $(MAKE) run
	cd 3_cuerpos && $(MAKE) run
clean:
	cd shocktube && $(MAKE) clean
	cd 3_cuerpos && $(MAKE) clean
