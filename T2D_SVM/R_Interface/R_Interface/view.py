from django.http import HttpResponse
from django.shortcuts import render,render_to_response
import uuid
import os
import json


def index(request):
    return render_to_response("T2D_svm.html")
def handle_upload_file(file,forlder_path):
    of = open(os.path.join(forlder_path,file.name), 'wb+')
    for chunk in file.chunks():
        of.write(chunk)
    of.close()

def upload_file(request):
    if request.method=="POST":
        #get file from request.FILES
        files=[]
        f_list=request.FILES
        for i in range(0,len(f_list)):
            files.append(f_list.get(str(i)))
        print f_list.get(str(0))
        print files
        if files == None:
            return HttpResponse(1)
        base_path=os.path.abspath(".")
        u_path = str(uuid.uuid1())
        folder_path = os.path.join(base_path,"upload/", u_path)
        if not os.path.exists(folder_path):
            os.mkdir(folder_path)
        for f in files:
            handle_upload_file(f,folder_path)
        return HttpResponse(str(folder_path))

def submit_job(request):
    if request.method=="POST":
        #get parameters from request.POST
        print request.POST
        #upload_folder_path=request.POST.get('folder_path',None)
        HbA1c = request.POST.get('HbA1c', None)
        Age = request.POST.get('Age', None)
        BMI = request.POST.get('BMI', None)
        TG = request.POST.get('TG',None)
        GLU = request.POST.get('GLU',None)
        GLU30 = request.POST.get('GLU30',None)
        GLU60 = request.POST.get('GLU60',None)
        GLU120 = request.POST.get('GLU120',None)
        CP = request.POST.get('CP',None)
        CP30 = request.POST.get('CP30',None)
        CP60 = request.POST.get('CP60',None)
        CP120 = request.POST.get('CP120',None)
        INS = request.POST.get('INS',None)
        INS30 = request.POST.get('INS30',None)
        INS60 = request.POST.get('INS60',None)
        INS120 = request.POST.get('INS120',None)

        # run r script by system call
        r_script_name="T2D_svm.html"
        base_path=os.path.abspath(".")
        #input1 = os.path.join(upload_folder_path,"miR_miR_corr.txt")
        #input2 = os.path.join(upload_folder_path,"mR_mR_corr.txt")
        #input3 = os.path.join(upload_folder_path,"miR_mR_corr.txt")
        #input4 = os.path.join(upload_folder_path,"miR_miR_corr_pvalue.txt")
        #input5 = os.path.join(upload_folder_path, "mR_mR_corr_pvalue.txt")
        #input6 = os.path.join(upload_folder_path, "miR_mR_corr_pvalue.txt")
        output_file_name="output.txt"
        output_path=os.path.join(base_path,"output","output")

        #the format of Command Line
        cmd="Rscript"+" "+r_script_name+" "+HbA1c+" "+Age+" "+BMI+" "+TG+" "+GLU+" "+GLU30+" "+GLU60+" "+GLU120+" "+CP +" "+CP30+" "+CP60+" "+CP120+" "+INS+" "+INS30+" "+INS30+" "+INS60+" "+INS120+" "+output_path
        code=os.system(cmd)
        print code
        if code==0:
            return_data=0
        else:
            return_data=1
        return HttpResponse(return_data)

def download(request):
    # download
    output_path = os.path.join(os.path.abspath("."),"output","output _ Table_module_list.txt")
    with open(output_path) as f:
        c = f.read()
    return HttpResponse(c)

def hello(request):
    return HttpResponse("hello")
