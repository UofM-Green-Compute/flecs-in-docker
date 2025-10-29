#!/bin/bash

( cat presentation.md ; echo '\newpage'; echo '# Full Source'; echo '```'; cat ../../src/asteroids_knn.cpp ; echo '```' ) | pandoc --pdf-engine=xelatex default-notes.yaml -  -o presentation-notes.pdf
