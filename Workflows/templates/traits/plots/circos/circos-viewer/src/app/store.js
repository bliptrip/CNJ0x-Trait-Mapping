import { configureStore } from '@reduxjs/toolkit';
import { loadState, saveState } from './localStorage'

import viewControllerReducer from '../features/viewController/viewControllerSlice';
import effectPlotReducer from '../features/effectPlot/effectPlotSlice';
import lodProfilePlotReducer from '../features/lodProfilePlot/lodProfilePlotSlice';

const persistedState = loadState();

export const store = configureStore({
  reducer: {
    viewController: viewControllerReducer,
    effectPlot: effectPlotReducer,
    lodProfilePlot: lodProfilePlotReducer
  },
  preloadedState: persistedState
});

store.subscribe(() => {
    saveState(store.getState());
});
