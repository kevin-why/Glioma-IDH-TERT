import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectPercentile, f_classif
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from tpot.builtins import ZeroCount
from tpot.export_utils import set_param_recursive

# NOTE: Make sure that the outcome column is labeled 'target' in the data file
tpot_data = pd.read_csv('/NAS/wg_zsx/notebooks/projects/wanghy/3sequence_IDH+TERT_features.csv', sep='COLUMN_SEPARATOR', dtype=np.float64)
features = tpot_data.drop('target', axis=1)
training_features, testing_features, training_target, testing_target = \
            train_test_split(features, tpot_data['target'], random_state=12)


exported_pipeline = make_pipeline(
    SelectPercentile(score_func=f_classif, percentile=96),
    ZeroCount(),
    MLPClassifier(alpha=0.0001, learning_rate_init=0.001)
)
# Fix random state for all the steps in exported pipeline
set_param_recursive(exported_pipeline.steps, 'random_state', 12)

exported_pipeline.fit(training_features, training_target)
results = exported_pipeline.predict(testing_features)
