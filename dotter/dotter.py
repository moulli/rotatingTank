from flask import Flask, render_template, url_for, request
import random
import os
import json
from collections import deque


app = Flask(__name__)


# Deque allowing to remember previous image
prev_img = deque(maxlen=2)

# Add labels to do if necessary
with open('labels_TODO.json', 'r+') as f:
	data = json.load(f)
	if len(data) == 0:
		image_files = os.listdir(os.path.join(os.getcwd(), 'static', 'dataset'))
		f.seek(0)
		json.dump(image_files, f, indent=4)
		f.truncate()


# Routes
@app.route('/')
@app.route('/index')
@app.route('/home')
def index():
	return render_template('index.html')

@app.route('/dotter', methods=['GET', 'POST'])
def dotter():
	# Load global variable prev_img and forms
	global prev_img
	# Deal with forms
	if request.method == 'POST':
		if 'previous' in request.form:
			if len(prev_img) == 2:
				prev_img.pop() # deleting unused image
			return render_template('dotter.html', img=prev_img[0], prev_len=len(prev_img))
		elif 'ptsval' in request.form:
			with open('labels.json', 'r+') as f:
				data = json.load(f)
				img_path = prev_img[-1]
				data[img_path[16:-4]] = request.form['ptsval']
				f.seek(0) # should reset file position to the beginning.
				json.dump(data, f, indent=4)
				f.truncate() # remove remaining part
			with open('labels_TODO.json', 'r+') as f:
				image_files = json.load(f)
				if img_path[16:] in image_files:
					image_files.remove(img_path[16:])
					f.seek(0)
					json.dump(image_files, f, indent=4)
					f.truncate()
			
	# Define new image
	with open('labels_TODO.json', 'r+') as f:
		image_files = json.load(f)
	image_file = url_for('static', filename='dataset/' + random.choice(image_files))
	prev_img.append(image_file)
	# Return dotter.html
	return render_template('dotter.html', img=image_file, prev_len=len(prev_img))


# Run app in debug mode
if __name__ == '__main__':
	app.run(debug=True, host='0.0.0.0', port=5000) 
	# host='0.0.0.0' is to make it accessible on the network through IP address
	# there is a bug whenever server is active on a port, then shut down then relaunched,
	# therefore, port=n allows to change port
	# NB: to connect locally, do not type 0.0.0.0:port, just 127.0.0.1:port
