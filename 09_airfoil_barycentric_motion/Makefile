.PHONY: all clean run

all:
	$(MAKE) -C backend

run: all
	./run.sh

clean:
	$(MAKE) -C backend clean
	rm -rf output
