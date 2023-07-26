rmarkdown.pandoc.to <- get0("rmarkdown.pandoc.to", ifnotfound = "docx")

is_word_output <- function () 
{
    return(grepl("docx", rmarkdown.pandoc.to))
}

is_html_output <- function () 
{
    return(grepl("html", rmarkdown.pandoc.to))
}

is_pdf_output <- function () 
{
    return(grepl("latex", rmarkdown.pandoc.to))
}