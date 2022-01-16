from flask import Flask, render_template, request, jsonify, send_file
import mysql.connector
from werkzeug.utils import secure_filename
import os
import re
from datetime import datetime
from subprocess import Popen
import time
import random
from config import *
 
app = Flask(__name__)
processes = []
         
app.secret_key = secret_key
emailregex = re.compile(email_regex_string)

   
def allowed_file(filename):
 return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def get_filename(filename):
    splitted = filename.lower().rsplit('.', 1)
    return getRandomUniqueString() + '.' + splitted[1]

def isValidEmail(email):
    return re.fullmatch(emailregex, email)

def getRandomUniqueString():
    return str(random.getrandbits(32)) + str(int(time.time()))
 
@app.route('/')
def index(): 
    return render_template('index.html')
 
@app.route("/upload",methods=["POST","GET"])
def upload():
    fileP = request.files['uploadPdbFile']
    filenameP = secure_filename(fileP.filename)
    fileX = request.files['uploadXtcFile']
    filenameX = secure_filename(fileX.filename)
    email = request.form['email']
    projectTitle = request.form['projectTitle']
    access = int(request.form['access'])
    resultId = getRandomUniqueString()

    if(not isValidEmail(email)):
        msg  = 'Invalid email address.'
        return jsonify({'htmlresponse': render_template('response.html', msg=msg, email=email)})

    if fileP and allowed_file(fileP.filename) and fileX and allowed_file(fileX.filename):
        cnx = mysql.connector.connect(user = 'root', password = 'admin12345678', host = '127.0.0.1', database = 'kavraki')
        filenamePathP = get_filename(filenameP)
        filenamePathX = get_filename(filenameX)
        pathP = os.path.join(upload_folder, filenamePathP)
        pathX = os.path.join(upload_folder, filenamePathX)
        fileP.save(pathP)
        fileX.save(pathX)
        today = datetime.today() 
        cur = cnx.cursor()
        sqlQuery = "INSERT INTO uploads (pdb_file_name, xtc_file_name,upload_time,email,title,resultId,access,username,password,pass_salt) VALUES (%s,%s,%s,%s,%s,%s,%s,'','','');"
        cur.execute(sqlQuery,[filenamePathP,filenamePathX,today,email,projectTitle,resultId,access])
        cmd = "python wfManager.py \"" + resultId + "\""
        cnx.commit()
        cur.close()
        cnx.close()
        print(cmd)
        file1 = open("templates/output_results/" + resultId + "_outputHTML.html","w")
        msg = f"Your analysis is being proccessed. Please refresh page later to check for the results or wait for a finished status notification that will be sent to {email}."
        file1.write(msg)
        file1.close()
        processes.append(Popen(cmd, shell=True))
        msg  = f'Files {fileP.filename} and {fileX.filename} are successfully uploaded to the database and in the process of analysis. Please check {email} for a notification when the analysis is finish. Your results will be available on the page: {domain_url}/view/{resultId}</a>'
    else:
        msg  = 'Invalid Upload - only pdb, xtc files are accepted.'
    return jsonify({'htmlresponse': render_template('response.html', msg=msg)})
 
@app.route('/view/<result_id>')
def view(result_id):
    return render_template('view.html', result_id=result_id)

@app.route('/getresult/<resultId>')
def getresult(resultId):
    msg  = 'Please wait.'
    path = "output_results/" + resultId + "_outputHTML.html"
    return jsonify({'htmlresponse': render_template(path, msg=msg)})

@app.route('/view/static/outputfiles/plots/<name>')
def get_image(name):
    path = "./static/outputfiles/plots/" + name
    return send_file(path, mimetype='image/gif')

if __name__ == "__main__":
    app.run(debug=True)