FILES=$(shell echo "src/*.md")
SUBFILES=$(shell echo "src/submodules/*.md")

make.jl: $(FILES) $(SUBFILES)
	julia make.jl
	formathtml .
	mkdir -p ../../temp
	cp -r build/. ../../temp/
	git checkout gh-pages
	rsync ../../temp/ ./
	git add .
	date=`date`
	git commit -m "updated docs ${date}"
	git checkout master
	rm -rf ../../temp

view: $(FILES) $(SUBFILES)
	julia make.jl
	firefox build/index.html

.PHONY: clean

clean: 
	rm -rf build/