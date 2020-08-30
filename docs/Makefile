FILES=$(shell echo "src/*.md")
SUBFILES=$(shell echo "src/submodules/*.md")

make.jl: $(FILES) $(SUBFILES)
	julia make.jl

.PHONY: clean

clean: 
	rm -rf build/