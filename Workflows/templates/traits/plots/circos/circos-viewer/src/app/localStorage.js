import DataFrame from 'dataframe-js';

export const loadState = () => {
    try {
        const serializedState = localStorage.getItem('circosState');
        if( serializedState === null ) {
            return undefined;
        }
        var savedState = JSON.parse(serializedState);
        return savedState;
    } catch(err) {
        return undefined;
    }
};

export const saveState = (state) => {
    try {
        const serializedState = JSON.stringify(state);
        localStorage.setItem('circosState', serializedState);
    } catch(err) {
        //Ignore error
    }
};

export const clearState = () => {
    try {
        localStorage.removeItem('circosState');
    } catch(err) {
        //Ignore error
    }
}
