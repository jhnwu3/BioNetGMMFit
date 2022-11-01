from flask import Flask, render_template, request, flash, redirect
from werkzeug.utils import secure_filename
from werkzeug.datastructures import  FileStorage
import subprocess
import os, shutil
import sys


def deleteFilesInDirectory(folder):
   for filename in os.listdir(folder):
      file_path = os.path.join(folder, filename)
      try:
         if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)
         elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
      except Exception as e:
         print('Failed to delete %s. Reason: %s' % (file_path, e))

app = Flask(__name__)

ALLFILES = []
INPUTS = {"config":"", "bngl":"", "timesteps":"", "truerates":"", "xData":"", "yData":[]}
OUTPUTS = {"estimates":"", "prediction":""}

# Get current path
path = os.getcwd()
# file Upload
UPLOAD_FOLDER = os.path.join(path, 'uploads')
OUTPUT_FOLDER = os.path.join(path, 'frontend/static')
# Make directory if uploads is not exists
if not os.path.isdir(UPLOAD_FOLDER):
    os.mkdir(UPLOAD_FOLDER)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Allowed extension you can set your own
ALLOWED_EXTENSIONS = set(['txt', 'pdf', 'png', 'jpg', 'jpeg', 'gif', 'csv', 'bngl'])
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/', methods=['GET'])
def index():
   return render_template('index.html', inputs = INPUTS, uploadedFiles = ALLFILES, outputs=OUTPUTS)

@app.route('/', methods=['POST'])
def upload_file():
    if request.method == 'POST':

        if 'files[]' not in request.files:
            flash('No file part')
            return redirect(request.url)

        files = request.files.getlist('files[]')

        for file in files:
            if file and allowed_file(file.filename):
               filename = secure_filename(file.filename)
               file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
               if filename not in ALLFILES:
                  ALLFILES.append(filename)
        flash('File(s) successfully uploaded')
        return redirect('/')

# configuration file
@app.route('/config/<file>', methods = ['GET'])
def select_config(file):
   if request.method == 'GET':
      INPUTS['config'] = file
      return redirect('/')
# model
@app.route('/model/<file>', methods = ['GET'])
def select_model(file):
   if request.method == 'GET':
      INPUTS['bngl'] = file
      return redirect('/')
# time steps
@app.route('/ts/<file>', methods = ['GET'])
def select_time(file):
   if request.method == 'GET':
      INPUTS['timesteps'] = file
      return redirect('/')
   
# x data
@app.route('/X/<file>', methods = ['GET'])
def select_initial(file):
   if request.method == 'GET':
      INPUTS['xData'] = file
      return redirect('/')

# true rate constants
# x data
@app.route('/tr/<file>', methods = ['GET'])
def select_rates(file):
   if request.method == 'GET':
      INPUTS['truerates'] = file
      return redirect('/')


# y data
@app.route('/Y/<file>', methods = ['GET'])
def select_evolved(file):
   if request.method == 'GET':
      if file not in INPUTS['yData']:
         INPUTS['yData'].append(file)
      return redirect('/')
   
@app.route('/reset', methods = ['GET'])
def reset():
   INPUTS['yData'] = []
   INPUTS['config'] = ""
   INPUTS['xData'] = ""
   INPUTS['timesteps'] = ""
   INPUTS['bngl'] = ""
   return redirect('/')

@app.route('/run', methods = ['POST'])
def run():
   if request.method == 'POST':
      allFull = True
      error = None
      for key, value in INPUTS.items():
         if key != 'truerates':
            if value == "" or value == []:
               allFull = True 
      
      deleteFilesInDirectory(UPLOAD_FOLDER + '/X')
      deleteFilesInDirectory(UPLOAD_FOLDER + '/Y')
      deleteFilesInDirectory(OUTPUT_FOLDER)
      shutil.move(UPLOAD_FOLDER + '/' + INPUTS['xData'], UPLOAD_FOLDER + '/X/' + INPUTS['xData'])
      for f in INPUTS['yData']:
         shutil.move(UPLOAD_FOLDER + '/' + f, UPLOAD_FOLDER + '/Y/' + f)
      
      if(allFull):
         if INPUTS['truerates'] == "":
            INPUTS['truerates'] = 'truerates.csv'
         subprocess.run(['./BNGMM', '-c', UPLOAD_FOLDER+ '/' + INPUTS['config'] ,
                         '-t', UPLOAD_FOLDER+ '/' + INPUTS['timesteps'],
                         '-m', UPLOAD_FOLDER+ '/' + INPUTS['bngl'], 
                         '-r', UPLOAD_FOLDER+ '/' + INPUTS['truerates'],
                         '-x', UPLOAD_FOLDER+'/X', 
                         '-y', UPLOAD_FOLDER+'/Y',
                         '-o', 'frontend/static/'], cwd=".")
         
         onlyfiles = [f for f in os.listdir(OUTPUT_FOLDER) if os.path.isfile(os.path.join(OUTPUT_FOLDER, f))]
         for f in onlyfiles:
            if "_estimates.png" in f: 
               OUTPUTS['estimates'] = 'static/' + f 
            if "_leastCostMoments.png" in f:
               OUTPUTS['prediction'] = 'static/' + f
         
         return redirect('/')
      else:
         error = "Missing All Required Inputs!" 
         flash(error)
         return redirect('/')
      
# def cygmm_output(path):
#    with open(path, 'r') as file:
#       data = file.readlines()
#       output = ""
#       for line in data:
#          output+= "<br>" + line + "</br>"
#       return output

# # come back to this..
# @app.route('/input', methods=['GET', 'POST'])
# def input():
#    print(request.method)
#    configExists = not inputs["config"] == ""
#    bnglExists = not inputs["bngl"] == ""
#    timestepsFileExists = not inputs["timesteps"] == ""
#    trueRatesFileExists = not inputs["truerates"] == ""
#    xDataExists = not inputs["xData"] == ""
#    yDataExists = not inputs["yData"] == []
#    if request.method == 'POST':
#       if request.form.get('Run') == 'Run':
#             # pass
#             print("Run")
#             subprocess.check_call(['./CyGMM'], cwd="../../", stdout=open("../../out.txt",'w'))
#             return render_template('output.html', output=cygmm_output("../../out.txt"))
#       else:
#             # pass # unknown
#             return render_template("index.html")
#    elif request.method == 'GET':
#       # return render_template("index.html")
#       print("No Post Back Call")
#    return render_template("index.html", config=configExists, bngl = bnglExists, ts = timestepsFileExists, tr = trueRatesFileExists, xData = xDataExists, yData = yDataExists)

# # @app.route('/upload')
# # def upload_file():
# #    return render_template('upload.html')
	
# # @app.route('/uploader/', methods = ['GET', 'POST'])
# # def upload_file():
# #    if request.method == 'POST':
# #       f = request.files['file']
# #       # print(request.form)
# #       # print(request.args)
# #       print(f.filename)
# #       f.save(secure_filename(f.filename))
# #       return 'file uploaded successfully'

# @app.route('/config', methods = ['GET', 'POST'])
# def upload_configuration():
#    if request.method == 'POST':
#       f = request.files['file']
#       # print(request.form)
#       # print(request.args)
#       print(f.filename)
#       if f.filename == '':
#          flash('Error no config file specified!')
#          return redirect("/")
      
#       inputs["config"] = f.filename
#       f.save(secure_filename(f.filename))
#       return render_template("config.html")
#    return render_template("config.html")

# def getEstimate(args, multi=False):
#    return args[args.index('-e') + 1]

# def getPrediction(args):
#    return args[args.index('-p') + 1]

if __name__ == '__main__':
   # outputPathing["estimates"] = "static/" + getEstimate(sys.argv)
   # outputPathing["prediction"] ="static/" + getPrediction(sys.argv)
   SECRET_KEY = os.urandom(24)
   app.secret_key = SECRET_KEY
   app.run(debug=True) 