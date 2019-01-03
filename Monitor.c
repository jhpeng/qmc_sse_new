#include <stdio.h>

void SetHTMLStyle(FILE* file)
{
    fprintf(file,"<style>\n");
    fprintf(file,"table {\n");
    fprintf(file,"\tfont-family: arial, sans-self;\n");
    fprintf(file,"\tborder-collapse: collapse;\n");
    fprintf(file,"\twidth: 100%%;\n");
    fprintf(file,"}\n");
    fprintf(file,"\n");
    fprintf(file,"td, th {\n");
    fprintf(file,"\tborder: 1px solid #dddddd;\n");
    fprintf(file,"\ttext-align: left;\n");
    fprintf(file,"\tpadding: 8px;\n");
    fprintf(file,"}\n");
    fprintf(file,"\n");
    fprintf(file,"tr:nth-child(even) {\n");
    fprintf(file,"\tbackground-color: #dddddd;\n");
    fprintf(file,"}\n");
    fprintf(file,"</style>\n");
}

void OutputHTML()
{
    FILE* file = fopen("test.html","w");
    fprintf(file,"<!DOCTYPE html>\n");
    fprintf(file,"<html>\n");
    fprintf(file,"<head>\n");
    SetHTMLStyle(file);
    fprintf(file,"</head>\n");
    fprintf(file,"<body>\n");
    fprintf(file,"\n");
    fprintf(file,"<h1>HTML Table</h1>\n");
    fprintf(file,"<h2>HTML Table</h2>\n");
    fprintf(file,"<h3>HTML Table</h3>\n");
    fprintf(file,"<table>\n");
    fprintf(file,"\t<tr>\n");
    fprintf(file,"\t\t<th>Company</th>\n");
    fprintf(file,"\t\t<th>Contact</th>\n");
    fprintf(file,"\t\t<th>Country</th>\n");
    fprintf(file,"\t</tr>\n");
    fprintf(file,"\t<tr>\n");
    fprintf(file,"\t\t<td>Company</td>\n");
    fprintf(file,"\t\t<td>Contact</td>\n");
    fprintf(file,"\t\t<td>Country</td>\n");
    fprintf(file,"\t</tr>\n");
    fprintf(file,"\t<tr>\n");
    fprintf(file,"\t\t<td>Company</td>\n");
    fprintf(file,"\t\t<td>Contact</td>\n");
    fprintf(file,"\t\t<td>Country</td>\n");
    fprintf(file,"\t</tr>\n");
    fprintf(file,"\t<tr>\n");
    fprintf(file,"\t\t<td>Company</td>\n");
    fprintf(file,"\t\t<td>Contact</td>\n");
    fprintf(file,"\t\t<td>Country</td>\n");
    fprintf(file,"\t</tr>\n");
    fprintf(file,"</table>\n");
    fprintf(file,"</body>\n");
    fprintf(file,"</html>\n");
}

#if 1
int main()
{
    OutputHTML();
}
#endif
