#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>

const float tau = 1.6180339887498949; // the golden ratio

typedef struct{
  double q[4];
}Quaternion;

int n_to_samples(int n){return 20*(n+5*pow(n,3));}

void quaternion_normalize(Quaternion *a)
{
  double abs = sqrt(pow(a->q[0],2) + pow(a->q[1],2) + pow(a->q[2],2) + pow(a->q[3],2));
  a->q[0] = a->q[0]/abs;
  a->q[1] = a->q[1]/abs;
  a->q[2] = a->q[2]/abs;
  a->q[3] = a->q[3]/abs;
}

double scalar_product_with_best_center(Quaternion quaternion, double * centers) {
  /* find closest center */
  double best_dist = 0.;
  int best_index = 0;
  double dist;
  for (int i = 0; i < 600; i++) {
    dist = ((quaternion.q[0] * centers[4*i+0] + quaternion.q[1] * centers[4*i+1] +
	     quaternion.q[2] * centers[4*i+2] + quaternion.q[3] * centers[4*i+3]) / 
	    sqrt(pow(quaternion.q[0], 2) + pow(quaternion.q[1], 2) + pow(quaternion.q[2], 2) + pow(quaternion.q[3], 2)));
    if (isinf(dist) || isnan(dist)) {
      printf("%d : %g %g %g %g\n", i, quaternion.q[0], quaternion.q[1], quaternion.q[2], quaternion.q[3]);
    }
    if (dist > best_dist) {
      best_dist = dist;
      best_index = i;
    }
  }

  /* calculate scalar product */
  double scalar_product = ((quaternion.q[0] * centers[4*best_index+0] + quaternion.q[1] * centers[4*best_index+1] +
			    quaternion.q[2] * centers[4*best_index+2] + quaternion.q[3] * centers[4*best_index+3]) /
			   sqrt(pow(quaternion.q[0], 2) + pow(quaternion.q[1], 2) + pow(quaternion.q[2], 2) + pow(quaternion.q[3], 2)));
  return scalar_product;
}

int generate_rotation_list(const int n, Quaternion *return_list, double *return_weights) {
  Quaternion *rotation_list = malloc(120*sizeof(Quaternion));

  for (int i = 0; i < 120; i++) {
    rotation_list[i].q[0] = 0.0;
    rotation_list[i].q[1] = 0.0;
    rotation_list[i].q[2] = 0.0;
    rotation_list[i].q[3] = 0.0;
  }

  /* first 16 */
  for (int i1 = 0; i1 < 2; i1++) {
    for (int i2 = 0; i2 < 2; i2++) {
      for (int i3 = 0; i3 < 2; i3++) {
	for (int i4 = 0; i4 < 2; i4++) {
	  rotation_list[8*i1+4*i2+2*i3+i4].q[0] = -0.5 + (double)i1;
	  rotation_list[8*i1+4*i2+2*i3+i4].q[1] = -0.5 + (double)i2;
	  rotation_list[8*i1+4*i2+2*i3+i4].q[2] = -0.5 + (double)i3;
	  rotation_list[8*i1+4*i2+2*i3+i4].q[3] = -0.5 + (double)i4;
	}
      }
    }
  }
  
  /* next 8 */
  for (int i = 0; i < 8; i++) {
    rotation_list[16+i].q[i/2] = -1.0 + 2.0*(double)(i%2);
  }

  /* last 96 */
  int it_list[12][4] = {{1,2,3,4},
		       {1,4,2,3},
		       {1,3,4,2},
		       {2,3,1,4},
		       {2,4,3,1},
		       {2,1,4,3},
		       {3,1,2,4},
		       {3,4,1,2},
		       {3,2,4,1},
		       {4,2,1,3},
		       {4,3,2,1},
		       {4,1,3,2}};

  

  for (int i = 0; i < 12; i++) {
    for (int j1 = 0; j1 < 2; j1++) {
      for (int j2 = 0; j2 < 2; j2++) {
	for (int j3 = 0; j3 < 2; j3++) {
	  rotation_list[24+8*i+4*j1+2*j2+j3].q[it_list[i][0]-1] = -0.5 + 1.0*(double)j1;
	  rotation_list[24+8*i+4*j1+2*j2+j3].q[it_list[i][1]-1] = tau*(-0.5 + 1.0*(double)j2);
	  rotation_list[24+8*i+4*j1+2*j2+j3].q[it_list[i][2]-1] = 1.0/tau*(-0.5 + 1.0*(double)j3);
	  rotation_list[24+8*i+4*j1+2*j2+j3].q[it_list[i][3]-1] = 0.0;
	}
      }
    }
  }
  
  /* get edges */
  /* all pairs of of vertices whose sum is longer than 3 is an edge */
  double dist2;
  int count = 0;
  double edge_cutoff = 3.24;

  int edges[720][2];
  for (int i = 0; i < 120; i++) {
    for (int j = 0; j < i; j++) {
      dist2 =
	pow(rotation_list[i].q[0] + rotation_list[j].q[0], 2) +
	pow(rotation_list[i].q[1] + rotation_list[j].q[1], 2) +
	pow(rotation_list[i].q[2] + rotation_list[j].q[2], 2) +
	pow(rotation_list[i].q[3] + rotation_list[j].q[3], 2);
      if (dist2 > edge_cutoff) {
	edges[count][0] = i;
	edges[count][1] = j;
	count++;
      }
    }
  }

  /* get faces */
  /* all pairs of edge and vertice whith a sum larger than 7.5 is a face */
  double face_cutoff = 7.5;
  int face_done[120];
  for (int i = 0; i < 120; i++) {face_done[i] = 0;}
  count = 0;
  int faces[1200][3];
  for (int i = 0; i < 720; i++) {
    face_done[edges[i][0]] = 1;
    face_done[edges[i][1]] = 1;
    for (int j = 0; j < 120; j++) {
      if (face_done[j]) {
	/* continue if the face has already been in a vertex,
	   including the current one */
	continue;
      }
      dist2 =
	pow(rotation_list[j].q[0] + rotation_list[edges[i][0]].q[0] +
	    rotation_list[edges[i][1]].q[0], 2) +
	pow(rotation_list[j].q[1] + rotation_list[edges[i][0]].q[1] +
	    rotation_list[edges[i][1]].q[1], 2) +
	pow(rotation_list[j].q[2] + rotation_list[edges[i][0]].q[2] +
	    rotation_list[edges[i][1]].q[2], 2) +
	pow(rotation_list[j].q[3] + rotation_list[edges[i][0]].q[3] +
	    rotation_list[edges[i][1]].q[3], 2);
      if (dist2 > face_cutoff) {
	faces[count][0] = edges[i][0];
	faces[count][1] = edges[i][1];
	faces[count][2] = j;
	count++;
      }
    }
  }

  /* get cells */
  /* all pairs of face and vertice with a sum larger than 13.5 is a cell */

  double cell_cutoff = 13.5;
  int cell_done[120];
  for (int i = 0; i < 120; i++) {cell_done[i] = 0;}
  count = 0;
  int cells[600][4];
  double cell_centers[4*600];
  double cell_center_norm;
  for (int j = 0; j < 120; j++) {
    cell_done[j] = 1;
    for (int i = 0; i < 1200; i++) {
      if (cell_done[faces[i][0]] || cell_done[faces[i][1]] || cell_done[faces[i][2]]) {
	continue;
      }
      dist2 =
	pow(rotation_list[faces[i][0]].q[0] + rotation_list[faces[i][1]].q[0] +
	    rotation_list[faces[i][2]].q[0] + rotation_list[j].q[0], 2) +
	pow(rotation_list[faces[i][0]].q[1] + rotation_list[faces[i][1]].q[1] +
	    rotation_list[faces[i][2]].q[1] + rotation_list[j].q[1], 2) +
	pow(rotation_list[faces[i][0]].q[2] + rotation_list[faces[i][1]].q[2] +
	    rotation_list[faces[i][2]].q[2] + rotation_list[j].q[2], 2) +
	pow(rotation_list[faces[i][0]].q[3] + rotation_list[faces[i][1]].q[3] +
	    rotation_list[faces[i][2]].q[3] + rotation_list[j].q[3], 2);
      if (dist2 > cell_cutoff) {
	cells[count][0] = faces[i][0];
	cells[count][1] = faces[i][1];
	cells[count][2] = faces[i][2];
	cells[count][3] = j;
	
	cell_centers[4*count+0] = (rotation_list[cells[count][0]].q[0] + rotation_list[cells[count][1]].q[0] +
				  rotation_list[cells[count][2]].q[0] + rotation_list[cells[count][3]].q[0]);
	cell_centers[4*count+1] = (rotation_list[cells[count][0]].q[1] + rotation_list[cells[count][1]].q[1] +
				  rotation_list[cells[count][2]].q[1] + rotation_list[cells[count][3]].q[1]);
	cell_centers[4*count+2] = (rotation_list[cells[count][0]].q[2] + rotation_list[cells[count][1]].q[2] +
				  rotation_list[cells[count][2]].q[2] + rotation_list[cells[count][3]].q[2]);
	cell_centers[4*count+3] = (rotation_list[cells[count][0]].q[3] + rotation_list[cells[count][1]].q[3] +
				  rotation_list[cells[count][2]].q[3] + rotation_list[cells[count][3]].q[3]);
	cell_center_norm = (sqrt(pow(cell_centers[4*count+0], 2) + pow(cell_centers[4*count+1], 2) +
				 pow(cell_centers[4*count+2], 2) + pow(cell_centers[4*count+3], 2)));
	cell_centers[4*count+0] /= cell_center_norm; cell_centers[4*count+1] /= cell_center_norm;
	cell_centers[4*count+2] /= cell_center_norm; cell_centers[4*count+3] /= cell_center_norm;
	count++;
      }
    }
  }

  /*variables used to calculate the weights */
  double alpha = acos(1.0/3.0);
  double f1 = 5.0*alpha/2.0/M_PI;
  double f0 = 20.0*(3.0*alpha-M_PI)/4.0/M_PI;
  double f2 = 1.0;
  double f3 = 1.0;

  int number_of_samples = n_to_samples(n);
  Quaternion *new_list = malloc(number_of_samples*sizeof(Quaternion));
  for (int i = 0; i < number_of_samples; i++) {
    new_list[i].q[0] = 0.;
    new_list[i].q[1] = 0.;
    new_list[i].q[2] = 0.;
    new_list[i].q[3] = 0.;
  }

  double *weights = malloc(number_of_samples*sizeof(double));
  double dist3;
  double scalar_product;

  /* copy vertices */
  for (int i = 0; i < 120; i++) {
    new_list[i].q[0] = rotation_list[i].q[0];
    new_list[i].q[1] = rotation_list[i].q[1];
    new_list[i].q[2] = rotation_list[i].q[2];
    new_list[i].q[3] = rotation_list[i].q[3];
    dist3 = pow(pow(new_list[i].q[0],2)+
		pow(new_list[i].q[1],2)+
		pow(new_list[i].q[2],2)+
		pow(new_list[i].q[3],2),(double)3/(double)2);
    scalar_product = scalar_product_with_best_center(new_list[i], cell_centers);
    weights[i] = f0*scalar_product/(double)number_of_samples/dist3;
  }

  /* split edges */
  int edges_base = 120;
  int edge_verts = (n-1);
  int index;
  for (int i = 0; i < 720; i++) {
    for (int j = 0; j < edge_verts; j++) {
      index = edges_base+edge_verts*i+j;
      for (int k = 0; k < 4; k++) {
	new_list[index].q[k] = 
	  (double)(j+1) / (double)(edge_verts+1) * rotation_list[edges[i][0]].q[k] +
	  (double)(edge_verts-j) / (double)(edge_verts+1) * rotation_list[edges[i][1]].q[k];
      }
      dist3 = pow(pow(new_list[index].q[0],2) + pow(new_list[index].q[1],2) +
		  pow(new_list[index].q[2],2) + pow(new_list[index].q[3],2), (double)3/(double)2);
      scalar_product = scalar_product_with_best_center(new_list[index], cell_centers);
      weights[index] = f1*scalar_product/(double)number_of_samples/dist3;
    }
  }

  /* split faces */
  int faces_base = 120 + 720*edge_verts;
  int face_verts = ((n-1)*(n-2))/2;
  double a,b,c;
  int kc;
  if (face_verts > 0) {
    for (int i = 0; i < 1200; i++) {
      count = 0;
      for (int ka = 2; ka < edge_verts+1; ka++) {
	for (int kb = 2; kb < edge_verts+1; kb++) {
	  if (ka + kb > edge_verts+1) {
	    kc = 2*(edge_verts+1)-ka-kb;
	    a = (double) (edge_verts + 1 - ka) / (double) (3*(edge_verts+1)-ka-kb-kc);
	    b = (double) (edge_verts + 1 - kb) / (double) (3*(edge_verts+1)-ka-kb-kc);
	    c = (double) (edge_verts + 1 - kc) / (double) (3*(edge_verts+1)-ka-kb-kc);
	    index = faces_base+face_verts*i+count;
	    for (int k = 0; k < 4; k++) {
	      new_list[index].q[k] =
		a * rotation_list[faces[i][0]].q[k] +
		b * rotation_list[faces[i][1]].q[k] +
		c * rotation_list[faces[i][2]].q[k];
	    }
	    dist3 = pow(pow(new_list[index].q[0],2) + pow(new_list[index].q[1],2) +
			pow(new_list[index].q[2],2) + pow(new_list[index].q[3],2), (double)3/(double)2);
	    scalar_product = scalar_product_with_best_center(new_list[index], cell_centers);
	    weights[index] = f2*scalar_product/(double)number_of_samples/dist3;
	    count++;
	  }
	}
      }
    }
  }

  /* split cells */
  int cell_base = 120 + 720*edge_verts + 1200*face_verts;
  int cell_verts = ((n-1)*(n-2)*(n-3))/6;
  double d;
  int kd;
  int debug_count = 0;
  if (cell_verts > 0) {
    for (int i = 0; i < 600; i++) { //600
      count = 0;
      for (int ka = 3; ka < edge_verts+1; ka++) {
	for (int kb = 3; kb < edge_verts+1; kb++) {
	  for (int kc = 3; kc < edge_verts+1; kc++) {
	    kd = 3*(edge_verts+1)-ka-kb-kc;
	    if (kd >= 3 && kd < edge_verts+1) {
	      a = (double) (edge_verts + 1 - ka) / (double) (4*(edge_verts+1)-ka-kb-kc-kd);
	      b = (double) (edge_verts + 1 - kb) / (double) (4*(edge_verts+1)-ka-kb-kc-kd);
	      c = (double) (edge_verts + 1 - kc) / (double) (4*(edge_verts+1)-ka-kb-kc-kd);
	      d = (double) (edge_verts + 1 - kd) / (double) (4*(edge_verts+1)-ka-kb-kc-kd);
	      index = cell_base+cell_verts*i+count;
	      for (int k = 0; k < 4; k++) {
		new_list[index].q[k] =
		  a*rotation_list[cells[i][0]].q[k] +
		  b*rotation_list[cells[i][1]].q[k] +
		  c*rotation_list[cells[i][2]].q[k] +
		  d*rotation_list[cells[i][3]].q[k];
	      }
	      dist3 = pow(pow(new_list[index].q[0],2) + pow(new_list[index].q[1],2) +
			  pow(new_list[index].q[2],2) + pow(new_list[index].q[3],2), (double)3/(double)2);
	      scalar_product = scalar_product_with_best_center(new_list[index], cell_centers);
	      weights[index] = f3*scalar_product/(double)number_of_samples/dist3;
	      count++;
	      debug_count++;
	    }
	  }
	}
      }
    }
  }

  free(rotation_list);

  /* prune list */
  /* Choose only quaternions with positive first element. If it is zero, the second element must be positive and so on. */
  const int number_of_rotations = n_to_samples(n);
  
  int counter = 0;
  
  Quaternion this_quaternion;
  int keep_this;
  for (int i = 0; i < number_of_rotations; ++i) {
    this_quaternion = new_list[i];
    keep_this = 0;
    if (this_quaternion.q[0] > 0) {
      keep_this = 1;
    } else if (this_quaternion.q[0] == 0.) {
      if (this_quaternion.q[1] > 0) {
	keep_this = 1;
      } else if(this_quaternion.q[1] == 0.) {
	if (this_quaternion.q[2] > 0) {
	  keep_this = 1;
	} else if(this_quaternion.q[2] == 0.) {
	  if (this_quaternion.q[3] > 0) {
	    keep_this = 1;
	  }
	}
      }
    }
    if (keep_this == 1) {
      return_list[counter] = this_quaternion;
      return_weights[counter] = weights[i];
      counter++;
    }
  }
  /* end prune */

  double weight_sum = 0.0;
  for (int i = 0; i < number_of_samples/2; i++) {
    weight_sum += return_weights[i];
  }
  for (int i = 0; i < number_of_samples/2; i++) {
    return_weights[i] /= weight_sum;
  }

  for (int i = 0; i < number_of_samples/2; i++) {
    quaternion_normalize(&return_list[i]);
  }
  
  return number_of_samples/2;
}

PyDoc_STRVAR(rotsampling__doc__, "rotsampling(sampling_n, return_weights)\n\nGenerate a uniform sampling of rotational space with a density of n (see [Loh2009]).");
static PyObject *rotsampling(PyObject *self, PyObject *args, PyObject *kwargs)
{
  int sampling_n;
  PyObject *return_weights_bool_obj = NULL;
  int return_weights_bool = 0;
  
  static char *kwlist[] = {"sampling_n", "return_weights", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i|O", kwlist, &sampling_n, &return_weights_bool_obj)) {
    return NULL;
  }

  if (sampling_n <= 0) {
    PyErr_SetString(PyExc_ValueError, "sampling_n must be positive.");
    return NULL;
  }

  if (return_weights_bool_obj == NULL) {
    return_weights_bool = 0;
  } else {
    if (PyObject_IsTrue(return_weights_bool_obj)) {
      return_weights_bool = 1;
    } else {
      return_weights_bool = 0;
    }
  }
  
  int number_of_rotations = n_to_samples(sampling_n)/2;
  int rotations_dim[] = {number_of_rotations, 4};
  PyObject *rotations_array = (PyObject *)PyArray_FromDims(2, rotations_dim, NPY_FLOAT64);
  double *rotations_raw = PyArray_DATA((PyArrayObject *)rotations_array);
  
  int weights_dim[] = {number_of_rotations};
  PyObject *weights_array = (PyObject *)PyArray_FromDims(1, weights_dim, NPY_FLOAT64);
  double *weights_raw = PyArray_DATA((PyArrayObject *)weights_array);

  generate_rotation_list(sampling_n, (Quaternion*)rotations_raw, weights_raw);

  if (return_weights_bool == 0) {
    Py_DECREF(weights_array);
    return rotations_array;
  } else {
    return Py_BuildValue("OO", rotations_array, weights_array);
  }
}

static PyMethodDef RotsamplingMethods[] = {
  {"rotsampling", (PyCFunction)rotsampling, METH_VARARGS|METH_KEYWORDS, rotsampling__doc__},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initrotsampling(void)
{
  import_array();
  PyObject *m = Py_InitModule3("rotsampling", RotsamplingMethods, "Sample rotation space");
  if (m == NULL)
    return;
}
