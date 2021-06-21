import { configureStore } from '@reduxjs/toolkit';
import viewControllerReducer from '../features/viewController/viewControllerSlice';

export const store = configureStore({
  reducer: {
    viewController: viewControllerReducer,
  },
});
