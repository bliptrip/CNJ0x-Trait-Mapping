all: render_html

gentextemplate:
		pandoc -D latex > template.tex
		sed -i -e 's#\\IfFileExists{footnote.sty}{\\usepackage{footnote}\\makesavenoteenv{longtable}}{}#% Commenting Out, as this breaks my longtable footnotes for bookdown -- \0#g' template.tex

render_office:
		@./renderOfficeDown.R && libreoffice _book/_main.docx
ro: render_office

render_word:
		@./renderWord.R && libreoffice _book/_main.docx
rw: render_word

package_word:
	rm -rf maule-ch1*
	mkdir -p maule-ch1/Data/publication
	cp _book/_main.docx maule-ch1/
	cp -r ../Data/publication maule-ch1/Data/
	zip -r maule-ch1.zip maule-ch1

pw: package_word


render_redoc:
		@./renderRedoc.R
rr: render_redoc

render_html:
		@./renderHTML.R && firefox _main.html
rh: render_html

render_pdf: gentextemplate
		@./renderPDF.R
rp: render_pdf

render_markdown:
		@./renderMD.R
rm: render_markdown

scpvh:
	scp -r -P 7262 _main.html maule@vaccinium.vcru.wisc.edu:/var/www/html/maule-ch1/index.html
	scp -r -P 7262 ../Data/publication maule@vaccinium.vcru.wisc.edu:/var/www/html/maule-ch1/Data/

