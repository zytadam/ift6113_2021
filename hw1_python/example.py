from trimesh import TriMesh
import meshplot as mp
import numpy as np
import os
import argparse
import math

def subdivision_method(inputTrimesh):
	# print("No subdivision implemented")
	oldVertexs = inputTrimesh.vs
	oldFaces = inputTrimesh.faces
	# print(oldFaces.dtype)

	tempFaces = np.c_[oldFaces[:,0], -np.ones((oldFaces.shape[0])),
					 oldFaces[:,1], -np.ones((oldFaces.shape[0])),
					 oldFaces[:,2], -np.ones((oldFaces.shape[0]))]
	
	# print(inputTrimesh.edges)
	# print(tempFaces)
	hes = inputTrimesh.get_halfedges()
	num = len(oldVertexs)

	edge_faces = []		#[ [ [v1, v2], faceIndex1, faceIndex2], ...]
	# Fill in temporary face matrix
	for (l,line) in enumerate(tempFaces):
		for (i,v) in enumerate(line):
			if v == -1:
				line[i] = num
				a = line[i-1]
				b = line[i+1 if i<5 else 0]
				ef = [[int(a), int(b)], l]

				hei = inputTrimesh.directed_edge2he_index((a, b))
				rhei = hes[hei].opposite_he
				if rhei != -1:
					iFace = hes[rhei].face
					pos = np.where(tempFaces[iFace] == b)[0][0] + 1
					tempFaces[iFace][pos] = num

					ef.append(iFace)
				num += 1
				edge_faces.append(ef)
	
	# print(edge_faces)

	# Build new face matrix
	newFaces = np.ndarray((0,3))
	newFaces = np.append(newFaces, np.c_[tempFaces[:,1:4]], axis = 0)
	newFaces = np.append(newFaces, np.c_[tempFaces[:,3:6]], axis = 0)
	newFaces = np.append(newFaces, np.c_[tempFaces[:,5], tempFaces[:,:2]], axis = 0)
	newFaces = np.append(newFaces, np.c_[tempFaces[:,1], tempFaces[:,3], tempFaces[:,5]], axis = 0)
	newFaces = newFaces.astype('int32')

	# print(tempFaces)
	# print(newFaces)

	# Relocate old vertexs
	newVertexs = np.ndarray((0,3))

	for (i, v) in enumerate(oldVertexs):
		neighboursIndex = inputTrimesh.vertex_vertex_neighbors(i)
		n = len(neighboursIndex)
		b = (5/8 - pow(3/8 + math.cos(2*math.pi/n)/4, 2))/n

		sum = np.zeros(3)
		for ineighbour in neighboursIndex:
			sum += oldVertexs[ineighbour]

		vp = (1 - n*b) * v + sum * b
		newVertexs = np.append(newVertexs, [vp], axis = 0)
	
	for vf in edge_faces:
		vp = np.zeros(3)

		if len(vf) == 3:
			oppos = np.setdiff1d(np.concatenate((oldFaces[vf[1]], oldFaces[vf[2]])), vf[0])
			for vi in vf[0]:
				vp += 3/8 * oldVertexs[vi]
			for vi in oppos:
				vp += 1/8 * oldVertexs[vi]

		elif len(vf) == 2:
			for vi in vf[0]:
				vp += 1/2 * oldVertexs[vi]

		newVertexs = np.append(newVertexs, [vp], axis = 0)
	# print(newVertexs)

	inputTrimesh.vs = newVertexs
	inputTrimesh.faces = newFaces
	inputTrimesh.topology_changed()

	return inputTrimesh

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="Run subdivision")
	parser.add_argument("--input", "-i", default="../input/cube.obj", help="path to input .obj")
	parser.add_argument("-n", default=3, type=int, help="number of iterations to perform")
	args = parser.parse_args()

	inputfile = args.input
	number_of_iterations = args.n

	mesh = TriMesh.FromOBJ_FileName(inputfile)

	os.makedirs("plots/", exist_ok=True)
	os.makedirs("output/", exist_ok=True)

	print("Saving a plot")
	mp.offline()
	p = mp.plot(mesh.vs, mesh.faces, c=mesh.vs[:,0], return_plot=True, filename='plots/test.html')

	for iteration in range(number_of_iterations):
		mesh = subdivision_method(mesh)

	mesh.write_OBJ("output/test_output.obj")
	print("Done")