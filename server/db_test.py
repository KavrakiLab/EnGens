import mysql.connector

cnx = mysql.connector.connect(user = 'root', password = 'admin12345678', host = '127.0.0.1', database = 'kavraki')
cnx.close()