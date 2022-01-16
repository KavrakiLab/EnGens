import nglview as nv
import matplotlib.pyplot as plt

class outputWriter():

    def __init__(self,resultId):
        self.resultId = resultId
        self.savePath = "./static/outputfiles/" + resultId + "/"
        self.HTML = ""
    
    def wf1intro(self):
        return
    
    def wf2intro(self):
        return
    
    def wf3intro(self):
        return

    def title(self, text):
        self.HTML += "<h1>" + str(text) + "</h1>"
        return 

    def write(self, text):
        self.HTML += "<p>" + str(text) + "<p>"
        return
    
    def writeFig(self, fig):
        #fig.savefig('path/to/save/image/to.png') 
        return
    
    def writeNGL(self,resultId, nglwidget, name):
        nglwidget.render_image()
        nglwidget._display_image()
        nv.write_html(self.savePath + name + ".html", nglwidget)
    
    def writePicture(self, text):
        self.HTML += "<img style='max-width: 100%'  src='" + str(text) + "' />"
        return

    def wf1step1(self):
        return
    
    def wf1step2(self):
        return
    
    def wf1step3(self):
        return

    def wf1step4(self):
        return

    def wf1step5(self):
        return

    def wf2step1(self):
        return
    
    def wf2step2(self):
        return
    
    def wf2step3(self):
        return

    def wf2step4(self):
        return

    def wf2step5(self):
        return

    def wf3step1(self):
        return
    
    def wf3step2(self):
        return
    
    def wf3step3(self):
        return

    def wf3step4(self):
        return
    
    def wf3step5(self):
        return