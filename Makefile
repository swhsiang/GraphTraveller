files = adj.c
project = proj3
map_file_5x5 = test/map5x5.txt
query_file_5x5 = test/query5x5.txt

compile: $(project)
	$(CC) -Wall -Werror -O3 $(files) -o $(project)

debug-compile: $(files)
	$(CC) -g $(files) -o $(project) 

llvm-compile: $(files)
	clang -g $(files) -o $(project)

test_sample:
	./$(project) $(map_file_5x5) $(query_file_5x5) 

clean:
	rm $(project) 
