all:
	make -C lib
	cp lib/cupyck.so cupyck/
	cp -r lib/parameters cupyck/

clean:
	make -C lib clean
	rm -f cupyck/cupyck.so
	rm -rf cupyck/parameters
