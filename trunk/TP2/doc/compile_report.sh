#!/bin/bash

rm ./report.log report.pdf report.toc report.dvi report.blg report.bbl report.aux;
pdflatex report.tex && pdflatex ./report.tex && ~/Adobe/Reader9/bin/acroread ./report.pdf
