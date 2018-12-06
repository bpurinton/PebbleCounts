# Install Pandoc on ubuntu
sudo apt install pandoc texlive-latex-recommended texlive-xetex texlive-luatex pandoc-citeproc etoolbox wkhtmltopdf

# Convert Markdown to PDF
(or any other format, just change filename extension)
d=$(date +%Y-%m-%d)
pandoc --number-sections --listings -H auto_linebreak_listings.tex \
    --variable papersize=a4paper --variable urlcolor=cyan \
    -s PebbleCounts_Manual_Nov2018_md_to_pdf.md -o PebbleCounts_Manual_Nov2018_md_to_pdf_${d}.pdf

cp PebbleCounts_Manual_Nov2018_md_to_pdf_${d}.pdf PebbleCounts_Manual.pdf
