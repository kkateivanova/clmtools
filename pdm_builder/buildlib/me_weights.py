import os, config, string
# these weights will weigh up some images where the eyes and mouth are in particular poses
# so that the SVM patches will be more generically discriminate

eyeweights = []
fi = open(os.path.join(config.data_folder, "training_hints/", "eyes_wide_open.csv"),"r")
for lines in fi:
	ewf = lines.strip().split(".")[:-1]
	eyeweights.append(string.join(ewf,".")+".bmp")
	eyeweights.append(string.join(ewf,".")+"_m.bmp")
fi.close()
fi = open(os.path.join(config.data_folder, "training_hints/", "eyes_closed.csv"),"r")
for lines in fi:
	ewf = lines.strip().split(".")[:-1]
	eyeweights.append(string.join(ewf,".")+".bmp")
	eyeweights.append(string.join(ewf,".")+"_m.bmp")
fi.close()
weights = dict.fromkeys([], eyeweights)

mouthweights = []
fi = open(os.path.join(config.data_folder, "training_hints/", "mouth.csv"),"r")
for lines in fi:
	ewf = lines.strip().split(".")[:-1]
	mouthweights.append(string.join(ewf,".")+".bmp")
	mouthweights.append(string.join(ewf,".")+"_m.bmp")
fi.close()
weights.update(dict.fromkeys([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17], mouthweights))
