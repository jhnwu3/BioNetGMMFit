from flask import Flask, render_template, request
import subprocess
import os

app = Flask(__name__)

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
   return render_template("index.html")

if __name__ == '__main__':
   app.run()