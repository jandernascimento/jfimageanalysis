#!/bin/bash

rm ./report.log report.pdf report.toc report.dvi report.blg report.bbl report.aux;
pdflatex report.tex && pdflatex ./report.tex && acroread -openInNewInstance ./report.pdf
