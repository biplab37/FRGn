FILES=$(shell echo "src/*.md")
SUBFILES=$(shell echo "src/submodules/*.md")

make.jl: $(FILES) $(SUBFILES)
	julia make.jl
	mkdir -p ../../temp
	cp -r build/ ../../temp/
	git checkout gh-pages
	cp -r ../../temp/ ./
	git add .
	git commit -m "update docs"
	git checkout master
	rm -rf ../../temp

view: $(FILES) $(SUBFILES)
	julia make.jl
	firefox build/index.html

.PHONY: clean

clean: 
	rm -rf build/