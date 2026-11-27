quarto<-Sys.which('quarto');if(!nzchar(quarto))quarto<-'C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe'
Sys.setenv(QUARTO_R=R.home('bin'),R_ENVIRON_USER=tempfile());setwd('report')
stopifnot(system2(quarto,c('render','Hospital_AMR_Report.qmd','--to','docx'))==0)
python<-Sys.getenv('HOSPITAL_AMR_PYTHON',Sys.which('python3'))
if(!nzchar(python))stop('Set HOSPITAL_AMR_PYTHON to a Python 3 executable for table formatting.')
stopifnot(system2(python,c('../scripts/style_report.py','Hospital_AMR_Report.docx'))==0)
