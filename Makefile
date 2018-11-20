files = adj.c
project = proj3

map_file_5x5 = test/map5x5.txt
query_file_5x5 = test/query5x5.txt
usa_data_set = test/usa.txt
usa_query_1 = test/usa1.txt
usa_query_10 = test/usa10.txt
usa_query_100 = test/usa100.txt

large_data_set = test/roadCA.txt
sample = test/shortestpath

compile: $(files)
	$(CC) -Wall -Werror -O3 $(files) -o $(project) -lm

compare: $(project) $(sample)
	@echo "Linux Only"
	@echo "time sample program"
	time ./$(sample) $(usa_data_set) $(usa_query_100)
	@echo "time the implementation"
	time ./$(project) $(usa_data_set) $(usa_query_100)


debug-compile: $(files)
	$(CC) -g $(files) -o $(project) 

llvm-compile: $(files)
	clang -g $(files) -o $(project)

test_sample:
	./$(project) $(map_file_5x5) $(query_file_5x5) 

clean:
	rm $(project) 
