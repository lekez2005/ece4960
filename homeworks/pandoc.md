pandoc --from gfm --to html  HW2.md | w3m -dump -T text/html -cols 200
pandoc --from gfm -t html5  HW2.md -o HW2.pdf --quiet
