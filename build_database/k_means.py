"""
Code for clustering sites into natural groups.
Based on https://datasciencelab.wordpress.com/2014/01/15/improved-seeding-for-clustering-with-k-means/
but with weighting (we act as if there are many points at the large wind/solar sites).
"""
from __future__ import division

import numpy as np
from collections import defaultdict

class KMeans(object):
    def __init__(self, K, X, size=None):
        if size is None:
            size = np.ones(len(X), dtype=float)
        self.X = np.array(X, dtype=float)
        # # choose up to K clusters, but not more than the number of unique elements in X
        # # this is an expensive test if x is big, so we don't do it
        # self.K = min(K, len(set(tuple(x) for x in self.X)))
        # choose K clusters
        self.K = K
        self.size = np.array(size, dtype=float)
        self.N = len(X)
        self.mu = np.array([])
        self.clusters = None
        self.method = None


    def _squared_distances(self):
        """
        return a matrix with one row for each entry in self.X and one column for
        each entry in self.mu, showing the squared distance between the corresponding
        entries in X and mu
        """
        D2 = np.array([
            np.sum((x - self.mu)**2, axis=1) for x in self.X
        ])
        # the next line is equivalent to above, but creates a 3D intermediate array,
        # which can be a problem for large datasets
        # (see http://scipy.github.io/old-wiki/pages/EricsBroadcastingDoc)
        # D2 = ((self.X[:,np.newaxis,:] - self.mu[np.newaxis,:,:])**2).sum(axis=2)
        return D2

    def _choose_next_center(self):
        # act as if there are a number of points at each location,
        # proportional to the specified sizes
        weights = np.copy(self.size)
        if len(self.mu) > 0:
            # one or more centers have already been selected;
            # give extra weight to the points far from the existing centers.
            weights *= self._squared_distances().min(axis=1)

        # choose a random point to add
        i = np.random.choice(self.X.shape[0], p=weights/np.sum(weights))
        x = self.X[i]
        if weights[0] is np.nan:
            raise ValueError('found NaN')

        if len(self.mu) > 0:
            # note: it's messy to keep re-creating mu as a numpy array,
            # but even if we used a list instead, that would get converted implicitly
            # to an array in the _squared_distances() calculation
            self.mu = np.append(self.mu, [x], axis=0)
        else:
            # no easy way to append a row array to an empty numpy array...
            self.mu = np.array([x])

    def init_centers(self, method='++'):
        self.method = method
        if method == 'random':
            # Initialize to K random centers (weighted by the project size at each location)
            self.mu = np.random.choice(self.X, size=self.K, p=self.size)
        else:   # method == '++'
            # initialize the centers using the k-means++ technique from
            # Arthur and Vassilvitskii (2007)
            # http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf
            while len(self.mu) < self.K:
                self._choose_next_center()

    # def _cluster_points(self):
    #     # self.clusters has one entry per cluster, which is a list of vectors
    #     # included in that cluster (the vectors themselves, not their IDs)
    #     self.clusters = defaultdict(list)
    #     best_mu_idx = self._squared_distances().argmin(axis=1)
    #     for i, x in enumerate(self.X):
    #         self.clusters[best_mu_idx[i]].append(x)
    #     self.cluster_id = best_mu_idx   # save cluster identifiers for plotting
    #
    # def _reevaluate_centers(self):
    #     for i in range(len(self.mu)):
    #         # self.clusters[i] is a list of rows of X that are in cluster i;
    #         # numpy calculates the mean across these rows as if this were an array
    #         self.mu[i] = np.mean(self.clusters[i], axis=0)

    def _group_mean(self, col):
        """
        col = self.X[:, 12]
        """
        return (
            np.bincount(self.cluster_id, weights=col*self.size, minlength=self.K)
            / np.bincount(self.cluster_id, weights=self.size, minlength=self.K)
        )

    def calculate_clusters(self):
        # find the id of the nearest cluster (0 to K-1) for each entry in X
        self.cluster_id = self._squared_distances().argmin(axis=1)
        # get weighted mean values of the X's that correspond to each cluster id
        self.mu = np.apply_along_axis(self._group_mean, axis=0, arr=self.X)

    def find_centers(self):
        while True:
            oldmu = np.copy(self.mu)
            # # Assign all points in X to clusters
            # self._cluster_points()
            # # Reevaluate centers
            # self._reevaluate_centers()
            # Reassign points in X to clusters and recalculate centers
            self.calculate_clusters()
            # check for convergence
            if np.array_equal(oldmu, self.mu):
                break
        return self.mu

    def plot(self):
        import matplotlib.pyplot as plt
        from scipy.spatial import Voronoi, voronoi_plot_2d
        # print "self.X:"
        # print self.X
        # print "self.size:"
        # print self.size
        fig = plt.figure(figsize=(10,8), dpi=96)
        ax = fig.add_subplot(111)
        # self.fig, self.ax = plt.subplots(1, 1)
        # self.fig.set_size_inches(10, 8)
        vor = Voronoi(self.mu)
        vor_plot = voronoi_plot_2d(vor, ax=ax)
        # remove the markers for each cluster and for the vertices
        ax.get_lines()[1].remove()
        #ax.get_lines()[0].remove()
        ax.scatter(x=self.X[:,0], y=self.X[:,1], c=self.cluster_id, s=self.size, alpha=0.75)
        ax.set_xlim(min(self.X[:,0]), max(self.X[:,0]))
        ax.set_ylim(min(self.X[:,1]), max(self.X[:,1]))
        # canvas.draw()   # may not be needed? see http://stackoverflow.com/questions/26783843/redrawing-a-plot-in-ipython-matplotlib

"""
Testing:

cell_capacities = np.load('cell_capacities.npy')
cell_cap_factors = np.load('cell_cap_factors.npy')
n_clusters = 4
self = KMeans(n_clusters, X=cell_cap_factors, size=cell_capacities)
self.init_centers()
self.find_centers()

"""
