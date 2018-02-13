#Makefile for building/packaging data & scripts for HTCondor submission

.PHONY: md


all: md

md: README.html

README.html: README.md
		pandoc -f markdown -t html -o $@ $<


install:
		R CMD BATCH install_packages.R
