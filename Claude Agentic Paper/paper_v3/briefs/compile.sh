#!/bin/bash
# Compile main.tex in /home/claude/work/paper_v3 ; creates placeholder PNGs for figures that do not exist yet.
cd /home/claude/work/paper_v3 || exit 1
for f in $(grep -v '^\s*%' main.tex | grep -o 'includegraphics\[[^]]*\]{[^}]*}' | sed 's/.*{\(.*\)}/\1/'); do
  base=$(basename "$f"); base=${base%.png}; base=${base%.pdf}
  if [ ! -f "figs/$base.png" ] && [ ! -f "figs/$base.pdf" ]; then
    python3 -c "
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig=plt.figure(figsize=(5.2,2.6)); fig.text(0.5,0.5,'PLACEHOLDER $base',ha='center',va='center',fontsize=14,color='red'); fig.savefig('figs/$base.png',dpi=72)"
    echo "placeholder created for $base"
  fi
done
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
bibtex main > /dev/null 2>&1
pdflatex -interaction=nonstopmode main.tex > /dev/null 2>&1
pdflatex -interaction=nonstopmode main.tex > build.log 2>&1
grep -A2 "^!" build.log | head -20
grep -i "undefined\|multiply" main.log | head
pdfinfo main.pdf | grep Pages
