import smtplib, ssl
import sys
from config import *

# sendEmail.py 
# receiver_email 
# result_id 
# projectTitle
# date 

password = "KavrakiLab123!"
port = 465
smtp_server = "smtp.gmail.com"
sender_email = "automated.message.from.the.lab@gmail.com"  # Enter your address

receiver_email =  sys.argv[1]
result_id = sys.argv[2]
projectTitle = sys.argv[3]
dateStamp = sys.argv[4]

message = f"""\
Subject: EnGens Result - Automated ensamble analysis results

Dear,
This is a notification to inform you that EnGens - Automated ensamble analysis for the {projectTitle} is finished. 
Results are available on the following webpage: {domain_url}/view/{result_id}.

This email is set on {dateStamp} for the result receiving.

Regards,
EnGens team
"""

context = ssl.create_default_context()
with smtplib.SMTP_SSL(smtp_server, port, context=context) as server:
    server.login(sender_email, password)
    server.sendmail(sender_email, receiver_email, message)