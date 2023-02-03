from flask import Flask, render_template, request, flash, redirect, copy_current_request_context, url_for, jsonify
from flask_executor import Executor
# from elasticsearch import Elasticsearch
# from flask_socketio import SocketIO, emit
# from celery import Celery
import subprocess
import json
import time
import random
from werkzeug.utils import secure_filename
from werkzeug.datastructures import  FileStorage
# import subprocess
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
app.config['EXECUTOR_TYPE'] = 'thread'
# app.config['EXECUTOR_MAX_WORKERS'] = 2
executor = Executor(app)


# with open('data/keys.json', 'r') as f:
#     keys = json.load(f)
#     ELASTIC_PASSWORD = keys['ELASTIC_PASSWORD']
#     CERT_FINGERPRINT = keys['CERT_FINGERPRINT']
    
# es = Elasticsearch(
#     "https://localhost:9200",
#     ssl_assert_fingerprint=(CERT_FINGERPRINT),
#     basic_auth=("elastic", ELASTIC_PASSWORD)
# )
# Celery configuration
# app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
# app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
# # Initialize Celery
# celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])

# socketio = SocketIO(app)
# socketio.init_app(app, cors_allowed_origins="*")
ALLFILES = []
INPUTS = {"config":"", "bngl":"","poi":"", "timesteps":"", "truerates":"", "xData":"", "yData":[]}
OUTPUTS = {"estimates":"", "prediction":""}
RUNNING = False
# Get current path
path = os.getcwd()
# file Upload
UPLOAD_FOLDER = os.path.join(path, 'uploads')
OUTPUT_FOLDER = os.path.join(path, 'frontend/static/')
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

# proteins of interest
@app.route('/poi/<file>', methods = ['GET'])
def select_proteins(file):
   if request.method == 'GET':
      INPUTS['poi'] = file
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
   INPUTS['truerates'] = ""
   INPUTS['poi'] = ""
   return redirect('/')

# @copy_current_request_context


@app.route('/run', methods = ['POST'])
def run():
   if request.method == 'POST':
      allFull = True
      error = None
      for key, value in INPUTS.items():
         if key != 'truerates' and key != 'poi':
            if value == "" or value == []:
               allFull = False 
               
      if(allFull):
         deleteFilesInDirectory(UPLOAD_FOLDER + '/X')
         deleteFilesInDirectory(UPLOAD_FOLDER + '/Y')
         deleteFilesInDirectory(OUTPUT_FOLDER)
         shutil.copyfile(UPLOAD_FOLDER + '/' + INPUTS['xData'], UPLOAD_FOLDER + '/X/' + INPUTS['xData'])
         for f in INPUTS['yData']:
            shutil.copyfile(UPLOAD_FOLDER + '/' + f, UPLOAD_FOLDER + '/Y/' + f)
         executor.submit_stored('bngmm', runBNGMM)
         # runBNGMM()
         # executor.add_default_done_callback(my_event)

         # if INPUTS['truerates'] == "":
         #    INPUTS['truerates'] = 'truerates.csv'
         # subprocess.run(['./BNGMM', '-c', UPLOAD_FOLDER+ '/' + INPUTS['config'] ,
         #                 '-t', UPLOAD_FOLDER+ '/' + INPUTS['timesteps'],
         #                 '-m', UPLOAD_FOLDER+ '/' + INPUTS['bngl'], 
         #                 '-r', UPLOAD_FOLDER+ '/' + INPUTS['truerates'],
         #                 '-x', UPLOAD_FOLDER+'/X', 
         #                 '-y', UPLOAD_FOLDER+'/Y',
         #                 '-o', 'frontend/static/'], cwd=".")
         # executor.submit(runBNGMM) 
      # socketio.emit('processing', json.dumps({'data': 'finished processing!'}))
      flash('Program Running! Please Refresh to See Results!')   
      # task = long_task.apply_async()
      return redirect('/result')
      # return jsonify({}), 202, {'Location': url_for('taskstatus',
      #                                             task_id=task.id)}
   else:
      error = "Missing All Required Inputs!" 
      flash(error)
      return redirect('/')


@app.route('/result')
def get_result():
   if not executor.futures.done('bngmm'):
      return render_template('result.html', done = False, outputs=OUTPUTS)
   future = executor.futures.pop('bngmm')
   return render_template('result.html', done = True, outputs=OUTPUTS)



# @executor.job
def runBNGMM():
   RUNNING = True
   if INPUTS['truerates'] == "":
      INPUTS['truerates'] = 'truerates.csv'
      
   bngmm = ['./BNGMM', '-c', UPLOAD_FOLDER + '/' + INPUTS['config'] ,
                     '-t', UPLOAD_FOLDER+ '/' + INPUTS['timesteps'],
                     '-m', UPLOAD_FOLDER+ '/' + INPUTS['bngl'], 
                     '-r', UPLOAD_FOLDER+ '/' + INPUTS['truerates'],
                     '-x', UPLOAD_FOLDER+'/X', 
                     '-y', UPLOAD_FOLDER+'/Y',
                     '-o', 'frontend/static/']
   
   ON_POSIX = 'posix' in sys.builtin_module_names
   # p = subprocess.Popen(bngmm, stdout=subprocess.PIPE, bufsize=1, close_fds=ON_POSIX, cwd= '.')
   p = subprocess.run(bngmm, cwd=".")
   # p.terminate()
   # executor.submit(my_event)
   # p.terminate()
   # p.wait()
   print(OUTPUT_FOLDER)
   onlyfiles = [f for f in os.listdir(OUTPUT_FOLDER) if os.path.isfile(os.path.join(OUTPUT_FOLDER, f))]
   print(onlyfiles)
   for f in onlyfiles:
      if "_estimates.png" in f: 
         OUTPUTS['estimates'] = '/static/'+ f 
      if "_leastCostMoments.png" in f:
         OUTPUTS['prediction'] = '/static/' + f
         
   print(OUTPUTS)
      # return render_template('index.html', inputs = INPUTS, uploadedFiles = ALLFILES, outputs=OUTPUTS)

if __name__ == '__main__':
   SECRET_KEY = os.urandom(24)
   app.secret_key = SECRET_KEY
   port = int(os.environ.get('PORT', 80))
   app.run(debug=True, host='0.0.0.0', port=port)
   # socketio.run(app, host='0.0.0.0', port=port, debug=True) 
   
   
   
   ### Defunct Code That I'm Keeping In Case, I ever want to revisit Celery
   
   # @celery.task(bind=True)
# def long_task(self):
#     """Background task that runs a long function with progress reports."""
#     verb = ['Starting up', 'Booting', 'Repairing', 'Loading', 'Checking']
#     adjective = ['master', 'radiant', 'silent', 'harmonic', 'fast']
#     noun = ['solar array', 'particle reshaper', 'cosmic ray', 'orbiter', 'bit']
#     message = ''
#     total = random.randint(10, 50)
#    #  for i in range(total):
#    #      if not message or random.random() < 0.25:
#    #          message = '{0} {1} {2}...'.format(random.choice(verb),
#    #                                            random.choice(adjective),
#    #                                            random.choice(noun))
#    #      self.update_state(state='PROGRESS',
#    #                        meta={'current': i, 'total': total,
#    #                              'status': message})
#    #      time.sleep(1)
#     runBNGMM()    
#     return {'current': 100, 'total': 100, 'status': 'Task completed!',
#             'result': 42}


# @app.route('/status/<task_id>')
# def taskstatus(task_id):
#     task = long_task.AsyncResult(task_id)
#     if task.state == 'PENDING':
#         response = {
#             'state': task.state,
#             'current': 0,
#             'total': 1,
#             'status': 'Pending...'
#         }
#     elif task.state != 'FAILURE':
#         response = {
#             'state': task.state,
#             'current': task.info.get('current', 0),
#             'total': task.info.get('total', 1),
#             'status': task.info.get('status', '')
#         }
#         if 'result' in task.info:
#             response['result'] = task.info['result']
#     else:
#         # something went wrong in the background job
#         response = {
#             'state': task.state,
#             'current': 1,
#             'total': 1,
#             'status': str(task.info),  # this is the exception raised
#         }
#     return jsonify(response)

# @socketio.on('long-running-event')
# def my_event(future):
#    global RUNNING
#    print("HEY IM HERE!")
#    # if executor.futures.done('run'):
#    #    executor.futures.pop('run')
#    #    onlyfiles = [f for f in os.listdir(OUTPUT_FOLDER) if os.path.isfile(os.path.join(OUTPUT_FOLDER, f))]
#    #    for f in onlyfiles:
#    #       if "_estimates.png" in f: 
#    #          OUTPUTS['estimates'] = 'static/' + f 
#    #       if "_leastCostMoments.png" in f:
#    #          OUTPUTS['prediction'] = 'static/' + f
#    #    socketio.emit('processing-finished', json.dumps({'data': 'finished processing!'}))
#    # else: 
#       socketio.emit('processing', json.dumps({'data': 'finished processing!'}))
      