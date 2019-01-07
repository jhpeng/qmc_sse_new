#include <stdio.h>
#include <string.h>

#include "DataStruct.h"
#include "Estimetor.h"

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

void OutputHTML(Observable* obs, SEPlaceHolder* placeholder, char* prefix)
{
    int isweep = placeholder->isweep+1;
    int nsweep = placeholder->nsweep;
    double ratio = (double)isweep/nsweep*100;
    char filename[128];
    sprintf(filename,"%s.html",prefix);
    FILE* file = fopen(filename,"w");
    fprintf(file,"<!DOCTYPE html>\n");
    fprintf(file,"<html>\n");
    fprintf(file,"<head>\n");
    SetHTMLStyle(file);
    fprintf(file,"</head>\n");
    fprintf(file,"<body>\n");
    fprintf(file,"\n");
    fprintf(file,"<h3>inverse temperature \t: %lf</h3>\n",placeholder->beta);
    fprintf(file,"<h3>length of sequence \t: %d</h3>\n",placeholder->ops->length);
    fprintf(file,"<h3>number of operator \t: %d</h3>\n",placeholder->ops->noo);
    fprintf(file,"<table>\n");
    fprintf(file,"\t<tr>\n");
    fprintf(file,"\t\t<th>obs name</th>\n");
    fprintf(file,"\t\t<th>mean</th>\n");
    fprintf(file,"\t\t<th>var</th>\n");
    fprintf(file,"\t\t<th>err</th>\n");
    fprintf(file,"\t</tr>\n");

    for(int i_obs=0;i_obs<obs->nobs;i_obs++){
        fprintf(file,"\t<tr>\n");
        fprintf(file,"\t\t<td>%s</td>\n",obs->obs_name[i_obs]);
        fprintf(file,"\t\t<td>%.4e</td>\n",obs->mean[i_obs]);
        fprintf(file,"\t\t<td>%.4e</td>\n",obs->var[i_obs]);
        fprintf(file,"\t\t<td>%.4e</td>\n",obs->err[i_obs]);
        fprintf(file,"\t</tr>\n");
    }
    fprintf(file,"</table>\n");
    fprintf(file,"<h3>|<----------------------%.2lf--------------------->| %d/%d</h3>\n",ratio,isweep,nsweep);
    fprintf(file,"</body>\n");
    fprintf(file,"</html>\n");
    fclose(file);
}

#if 0
int main()
{
    OutputHTML();
}
#endif
