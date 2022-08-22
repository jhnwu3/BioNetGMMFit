from flask import Flask, render_template, request, flash, redirect
from werkzeug.utils import secure_filename
from werkzeug.datastructures import  FileStorage
import subprocess
import os

app = Flask(__name__)
inputs = {"config":"", "bngl":"", "timesteps":"", "truerates":"", "xData":"", "yData":[]}


def cygmm_output(path):
   with open(path, 'r') as file:
      data = file.readlines()
      output = ""
      for line in data:
         output+= "<br>" + line + "</br>"
      return output

@app.route('/', methods=['GET', 'POST'])
def index():
   print(request.method)
   configExists = not inputs["config"] == ""
   bnglExists = not inputs["bngl"] == ""
   timestepsFileExists = not inputs["timesteps"] == ""
   trueRatesFileExists = not inputs["truerates"] == ""
   xDataExists = not inputs["xData"] == ""
   yDataExists = not inputs["yData"] == []
   if request.method == 'POST':
      if request.form.get('Run') == 'Run':
            # pass
            print("Run")
            subprocess.check_call(['./CyGMM'], cwd="../../", stdout=open("../../out.txt",'w'))
            return render_template('output.html', output=cygmm_output("../../out.txt"))
      else:
            # pass # unknown
            return render_template("index.html")
   elif request.method == 'GET':
      # return render_template("index.html")
      print("No Post Back Call")
   return render_template("index.html", config=configExists, bngl = bnglExists, ts = timestepsFileExists, tr = trueRatesFileExists, xData = xDataExists, yData = yDataExists)

# @app.route('/upload')
# def upload_file():
#    return render_template('upload.html')
	
# @app.route('/uploader/', methods = ['GET', 'POST'])
# def upload_file():
#    if request.method == 'POST':
#       f = request.files['file']
#       # print(request.form)
#       # print(request.args)
#       print(f.filename)
#       f.save(secure_filename(f.filename))
#       return 'file uploaded successfully'

@app.route('/config', methods = ['GET', 'POST'])
def upload_configuration():
   if request.method == 'POST':
      f = request.files['file']
      # print(request.form)
      # print(request.args)
      print(f.filename)
      if f.filename == '':
         flash('Error no config file specified!')
         return redirect("/")
      
      inputs["config"] = f.filename
      f.save(secure_filename(f.filename))
      return render_template("config.html")
   return render_template("config.html")


if __name__ == '__main__':
   SECRET_KEY = os.urandom(24)
   app.secret_key = SECRET_KEY
   app.run(debug=True) 