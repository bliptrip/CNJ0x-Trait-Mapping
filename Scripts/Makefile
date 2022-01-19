all: render_html

render_office:
		@./renderOfficeDown.R && libreoffice _book/_main.docx

ro: render_office

render_word:
		@./renderWord.R && libreoffice _book/_main.docx
rw: render_word

render_html:
		@./renderHTML.R && firefox _main.html
rh: render_html

render_pdf:
		@./renderPDF.R
rp: render_pdf

scpvh:
	scp -r -P 7262 _main.html maule@vaccinium.vcru.wisc.edu:/var/www/html/maule-ch1/index.html
	scp -r -P 7262 ../Data/publication maule@vaccinium.vcru.wisc.edu:/var/www/html/maule-ch1/Data/

