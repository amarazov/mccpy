FORT=gfortran
F2PY=f2py

all: ssobol.so sfaur.so

#main: main.o libssobol.so libsfaur.so
#	$(FORT) -o main main.o libssobol.so
	
#main.o: main.f
#	$(FORT) -c main.f
	
#libsfaur.so: sfaur.o
#	$(FORT) -shared -o libsfaur.so sfaur.o

#sfaur.o: sfaur.f
#	$(FORT)  -fno-align-commons -c -fPIC sfaur.f
	
sfaur.so:  sfaur.pyf sfaur.f sfaur.so
	$(F2PY) -c -m sfaur sfaur.pyf sfaur.f
	mv sfaur.so ../Python/
	
#libssobol.so: ssobol.o
#	$(FORT) -shared -o libssobol.so ssobol.o

#ssobol.o: ssobol.f
#	$(FORT)  -fno-align-commons -c -fPIC ssobol.f
	
ssobol.so:  ssobol.pyf ssobol.f
	$(F2PY) -c -m ssobol ssobol.pyf ssobol.f
	mv ssobol.so ../Python/

clean:
	rm -rf *o main
