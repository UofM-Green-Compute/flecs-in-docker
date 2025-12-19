
all: presentation.pdf

presentation.pdf : presentation.md default.yaml
	pandoc -t beamer default.yaml  presentation.md -o presentation.pdf
	./build-notes.sh

clean:
	rm -f *pdf

