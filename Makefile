CFLAGS=-Wall -O3

eQCM: bin/eQCM

bin/eQCM: .obj/Itemset.o .obj/Transactionset.o .obj/CartesianProduct.o .obj/WGCStatic.o .obj/WG.o .obj/CartesianProductDb.o .obj/BitDb.o .obj/main.o
	g++ $(CFLAGS) $^ -o bin/eQCM

.obj/%.o: src/%.cpp
	g++ $(CFLAGS) $< -c -o $@

clean:
	rm -f .obj/*.o 
