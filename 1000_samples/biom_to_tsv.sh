#!/bin/bash

for FILE in *.tsv; do sed -i '1d' $FILE ; sed -i  '1s/^..//' $FILE ;  done
