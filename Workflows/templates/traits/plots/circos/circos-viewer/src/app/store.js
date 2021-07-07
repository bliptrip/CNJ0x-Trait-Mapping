import { configureStore } from '@reduxjs/toolkit';
import { loadState, saveState } from './localStorage'

import viewReducer from '../features/view/viewSlice';
import viewControllerReducer from '../features/viewController/viewControllerSlice';
import effectPlotReducer from '../features/effectPlot/effectPlotSlice';
import corrPlotReducer from '../features/corrPlot/corrPlotSlice';
import lodProfilePlotReducer from '../features/lodProfilePlot/lodProfilePlotSlice';
import blupTableSliceReducer from '../features/blupTable/blupTableSlice';
import blupTableGridSliceReducer from '../features/blupTableGrid/blupTableGridSlice';

const persistedState = loadState();

export const store = configureStore({
  reducer: {
    view: viewReducer,
    viewController: viewControllerReducer,
    effectPlot: effectPlotReducer,
    corrPlot: corrPlotReducer,
    lodProfilePlot: lodProfilePlotReducer,
    blupTable: blupTableSliceReducer,
    blupTableGrid: blupTableGridSliceReducer
  },
  preloadedState: persistedState
});

store.subscribe(() => {
    saveState(store.getState());
});
