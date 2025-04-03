import numpy as np
import csv
import os
import re


def ReadCFXExport(fileName, hasFaces=False):
	if os.path.isfile(fileName):
		headerRows = 6
		with open(fileName) as dataFile:
			for i in range(headerRows-1): dataFile.readline()
			reader = csv.reader(dataFile)
			
			headers = []
			isVector = []
			for header in next(reader):
				h = re.sub(r'\[[^>]*\]', '', header).lstrip().rstrip()
				headers.append(h)
				
				isVector.append(('Gradient' in h) or ('Velocity' == h) or ('Velocity.Trnavg' == h) or ('Normal' in h))
			
			isFaces = False
			dataRowCount = 0
			for line in reader:
				if len(line) == 0:
					isFaces = True
					break
				dataRowCount += 1
		
		formats = np.where(isVector,'3f4','f4').tolist()
		if 'Node Number' in headers:
			formats[headers.index('Node Number')] = 'u4'
		if 'Velocity.Divergence' in headers:
			formats[headers.index('Velocity.Divergence')] = 'f4'
		
		'''
		vgrad = 'Velocity.Gradient'
		if 'Velocity u.Gradient' in headers and 'Velocity v.Gradient' in headers and 'Velocity w.Gradient' in headers:
			formats[headers.index('Velocity u.Gradient')] = '9f4'
			del(formats[headers.index('Velocity v.Gradient')])
			del(formats[headers.index('Velocity w.Gradient')])
			
			headers[headers.index('Velocity u.Gradient')] = vgrad
			del(headers[headers.index('Velocity v.Gradient')])
			del(headers[headers.index('Velocity w.Gradient')])
		
		print(len(headers))
		print(len(formats))
		for h, f in zip(headers,formats):
			print('{}	{}'.format(h,f))
		'''
		
		# print(headers, len(headers))
		# print(formats, len(formats))
		
		dtypes = {'names' : headers, 'formats' : formats}
		data = np.genfromtxt(fileName,
							delimiter=',',
							skip_header=headerRows,
							max_rows=dataRowCount,
							dtype=dtypes)
		
		'''
		if vgrad in headers:
			data[vgrad] = data[vgrad].reshape(-1,3,3)

		
		faces = np.genfromtxt(fileName,
					delimiter=',',
					skip_header=headerRows+dataRowCount+2,dtype='u4')
		'''
		
		return data
	else:
		print('No file: {}'.format(fileName))