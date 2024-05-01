import numpy as np
import random
from sklearn.tree import DecisionTreeClassifier
from sklearn.base import BaseEstimator
import multiprocessing as mp


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

class RandomForestClassifierCustom(BaseEstimator):
    def __init__(self, n_estimators=10, max_depth=None, max_features=None, random_state=SEED):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []
        self.random_states_by_tree = []

    def _fit_tree(self, X, y, random_state, feat_ids):
        np_random = np.random.RandomState(random_state)
        indices = np_random.choice(range(X.shape[0]), size=X.shape[0], replace=True)
        X_bootstrap = X[indices][:, feat_ids]
        y_bootstrap = y[indices]

        tree = DecisionTreeClassifier(max_depth=self.max_depth, max_features=self.max_features, random_state=random_state)
        tree.fit(X_bootstrap, y_bootstrap)
        return tree

    def fit(self, X, y, n_jobs=None):
        self.classes_ = sorted(np.unique(y))

        if n_jobs is None:
            n_jobs = mp.cpu_count()
        n_jobs = min(n_jobs, mp.cpu_count())

        args = [(X, y, self.random_state + i, np.random.choice(range(X.shape[1]), size=self.max_features, replace=False)) for i in range(self.n_estimators)]
        
        with mp.Pool(processes=n_jobs) as pool:
            self.trees = pool.starmap(self._fit_tree, args)
            self.feat_ids_by_tree = [arg[3] for arg in args]
            self.random_states_by_tree = [arg[2] for arg in args]

        return self

    def _predict_proba_tree(self, X, tree, feat_ids):
        return tree.predict_proba(X[:, feat_ids])

    def predict_proba(self, X, n_jobs=None):
        if n_jobs is None:
            n_jobs = mp.cpu_count()
        n_jobs = min(n_jobs, mp.cpu_count())

        with mp.Pool(processes=n_jobs) as pool:
            args = [(X, tree, self.feat_ids_by_tree[i]) for i, tree in enumerate(self.trees)]
            all_probas = pool.starmap(self._predict_proba_tree, args)

        probas = np.sum(all_probas, axis=0) / len(self.trees)
        return probas

    def predict(self, X, n_jobs=None):
        probas = self.predict_proba(X, n_jobs=n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions
