# Install Pandoc on ubuntu
sudo apt install pandoc texlive-latex-recommended texlive-xetex texlive-luatex pandoc-citeproc etoolbox wkhtmltopdf

# Convert Markdown to PDF
(or any other format, just change filename extension)
d=$(date +%Y-%m-%d)
pandoc --number-sections --listings -H auto_linebreak_listings.tex \
    --variable papersize=a4paper --variable urlcolor=cyan \
    -s From_PC_to_GMT_Maps.md -o From_PC_to_GMT_Maps_${d}.pdf

cp From_PC_to_GMT_Maps_${d}.pdf From_PC_to_GMT_Maps.pdf
