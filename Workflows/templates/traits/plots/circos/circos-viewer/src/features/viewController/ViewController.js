import React, { useEffect, useRef, useState, useCallback } from 'react';
import { useSelector, shallowEqual, useDispatch } from 'react-redux';
import { useDrag, useDrop } from 'react-dnd';
import {selectList,
        selectListItem,
        selectDisplayInteractions,
        selectDisplayTrackLabels,
        selectQTLModelCount,
        selectQTLConsensus,
        selectQTLMethod,
        setQTLModelCount,
        setQTLConsensus,
        setQTLMethod,
        setList,
        setListItemChecked,
        setDisplayInteractions,
        setDisplayTrackLabels,
        MAX_QTL_MODEL_COUNT
} from './viewControllerSlice';
import {clearState} from '../../app/localStorage';
//import styles from './ViewController.module.css';
import update from 'immutability-helper';

/* Material-UI Core and Components */
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Checkbox from '@material-ui/core/Checkbox';
import Divider from '@material-ui/core/Divider';
import FormGroup from '@material-ui/core/FormGroup';
import Grid from '@material-ui/core/Grid';
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormControl from '@material-ui/core/FormControl';
import Input from '@material-ui/core/Input';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import Modal from '@material-ui/core/Modal';
import Slider from '@material-ui/core/Slider';
import Switch from '@material-ui/core/Switch';
import Tooltip from '@material-ui/core/Tooltip';
import Typography from '@material-ui/core/Typography';

import { makeStyles } from '@material-ui/core/styles';

/* Material-UI Icons */
import AllInclusiveIcon from '@material-ui/icons/AllInclusive';
import ClearIcon from '@material-ui/icons/Clear';
import PlayForWorkIcon from '@material-ui/icons/PlayForWork';
import RotateLeftIcon from '@material-ui/icons/RotateLeft';

const ItemTypes = {
    VIEW_ELEMENT: 'ViewElement',
};

const useStyles = makeStyles((theme) => ({
    input: {
        width: 50,
    },
    button: {
        margin: theme.spacing(0.5),
    },
    paper: {
        position: 'absolute',
        top: "50%",
        left: "50%",
        transform: 'translate("-50%", "-50%")',
        backgroundColor: theme.palette.background.paper,
        border: '2px solid #000',
        boxShadow: theme.shadows[5],
        padding: theme.spacing(2, 4, 3),
    },
}));

const ConfigurationResetor = () => {
    const classes             = useStyles();
    const [open, setOpen]     = useState(false);

    const handleReset = (event) => {
        clearState();
        setOpen(true);
    }

    const handleClose = (event) => {
        setOpen(false);
    }

    return (
        <React.Fragment>
            <ListItem>
                <Tooltip title="Reset Configuration State to Defaults" arrow>
                    <Button
                        variant="contained"
                        color="primary"
                        className={classes.button}
                        startIcon={<RotateLeftIcon />}
                        onClick={handleReset}
                    >
                            Reset State
                    
                    </Button>
                </Tooltip>
                <Modal
                    open={open}
                    onClose={handleClose}
                    aria-labelledby="modal-reset"
                    aria-describedby="modal-description"
                >
                    <div className={classes.paper}>
                        Refresh browser to view changes.
                    </div>
                </Modal>
            </ListItem>
        </React.Fragment>
    );
};

const QTLModelCountSelector = () => {
    const classes             = useStyles();
    const dispatch            = useDispatch();
    const qtlCount            = useSelector(selectQTLModelCount);
    const [value, setValue]   = useState(qtlCount);

    const handleSliderChange = (event, newValue) => {
        setValue(newValue);
    };

    const handleInputChange = (event) => {
        setValue(event.target.value === '' ? '' : Number(event.target.value));
    };

    const handleSubmit = (event) => {
        dispatch(setQTLModelCount(value));
    }

    const handleBlur = () => {
        if (value < 1) {
            setValue(0);
        } else if (value > MAX_QTL_MODEL_COUNT) {
            setValue(MAX_QTL_MODEL_COUNT);
        }
    };

    return (
        <React.Fragment>
            <ListItem>
                <Tooltip title="Maximum number of QTLs to show per model-trait, sorted by percent variance explained." arrow>
                    <Typography id="qtl-count-slider" gutterBottom>
                        Max Number QTLs
                    </Typography>
                </Tooltip>
            </ListItem>
            <ListItem>
                <Grid container spacing={1} justifyContent="space-between" alignItems="center">
                    <Grid item>
                        <Tooltip title="Submit QTL Model Count" arrow>
                            <Button
                                variant="contained"
                                color="primary"
                                className={classes.button}
                                startIcon={<PlayForWorkIcon />}
                                onClick={handleSubmit}
                            />
                        </Tooltip>
                    </Grid>
                    <Grid item xs>
                        <Slider
                            value={typeof value === 'number' ? value : MAX_QTL_MODEL_COUNT}
                            onChange={handleSliderChange}
                            aria-labelledby="qtl-count-slider"
                            valueLabelDisplay="auto"
                            step={1}
                            marks
                            min={1}
                            max={MAX_QTL_MODEL_COUNT}
                        />
                    </Grid>
                    <Grid item>
                        <Input
                        className={classes.input}
                        value={value}
                        margin="dense"
                        onChange={handleInputChange}
                        onBlur={handleBlur}
                        inputProps={{
                            step: 1,
                            min: 1,
                            max: MAX_QTL_MODEL_COUNT,
                            type: 'number',
                            'aria-labelledby': 'qtl-count-slider',
                        }}
                        />
                    </Grid>
                </Grid>
            </ListItem>
        </React.Fragment>
    );
};


const QTLMethodSelector = () => {
    const dispatch      = useDispatch();
    const qtlMethod     = useSelector(selectQTLMethod);

    const handleChange  = (event) => {
        dispatch(setQTLMethod(event.target.value));
    };

    return (
        <List>
            <ListItem>
                <ListItemText primary="QTL Method" />
            </ListItem>
            <ListItem>
                <FormControl component="fieldset">
                    <RadioGroup row aria-label="qtl-method-selector" name="qtl-method-selector" defaultValue={qtlMethod} onChange={handleChange}>
                    <FormControlLabel
                        value="scanone"
                        control={<Radio color="primary" />}
                        label="r/QTL scanone()"
                        labelPlacement="end"
                    />
                    <FormControlLabel
                        value="stepwiseqtl"
                        control={<Radio color="primary" />}
                        label="r/QTL stepwiseqtl()"
                        labelPlacement="end"
                    />
                    </RadioGroup>
                </FormControl>
            </ListItem>
        </List>
    );
};

const QTLConsensusSwitch = () => {
    const dispatch      = useDispatch();
    const qtlConsensus  = useSelector(selectQTLConsensus);

    const handleChange  = (event) => {
        dispatch(setQTLConsensus(event.target.checked));
    };

    return( <ListItem>
                <FormGroup>
                    <FormControlLabel
                        control={
                                    <Switch
                                        checked={qtlConsensus}
                                        onChange={handleChange}
                                        color="primary"
                                    />
                                }
                        label="Consensus"
                    />
                </FormGroup>
            </ListItem> );
};

const ListViewController = ({listName, listType}) => {
    const classes                                = useStyles();
    const dispatch                               = useDispatch();
    const dispatchFn                             = setList(listType);
    const listItems                              = useSelector(selectList(listType), shallowEqual);
    const [listItemsState,setListItemsState]     = useState(JSON.parse(JSON.stringify(listItems)));
    const [globalSelect, setGlobalSelect]        = useState(false);

    useEffect(() => {
        setListItemsState(listItems);
    }, [listItems]);


    const handleListViewSubmit = (event) => {
        dispatch(dispatchFn(listItemsState)); 
    };

    const setListElementChecked = (id,checked) => {
        setListItemsState(
            listItemsState.map( e => {
                if( e.id == id ) {
                    return( { ...e, enabled: checked} );
                } else {
                    return e; //Unchanged
                }
            })
        );
    };

    const moveListViewElement = useCallback((dragIndex, hoverIndex) => {
        const item = listItemsState[dragIndex];
        setListItemsState(update(listItemsState, {
            $splice: [
                [dragIndex, 1],
                [hoverIndex, 0, item],
            ],
        }));
    }, [listItemsState, dispatch, dispatchFn]);

    const selectAll = (event) => {
        setListItemsState(listItemsState.map( e => ({...e, enabled: true}) ));
    }

    const selectNone = (event) => {
        setListItemsState(listItemsState.map( e => ({...e, enabled: false}) ));
    }

    const renderGlobalSelect = () => {
        const listItemsLength = listItemsState.length;

        if(listItemsLength > 3) { //Only render 'select all' and 'select none' if above a given threshold of entries in list
            const listItemsEnabledLength = listItemsState.filter( e => (e.enabled) ).length;
            const disableSelectAll = (listItemsEnabledLength === listItemsLength);
            const disableSelectNone = (listItemsEnabledLength === 0);
            return(
                <div>
                    <Tooltip title="Select All" arrow>
                        <Button
                            variant="contained"
                            color="primary"
                            className={classes.button}
                            startIcon={<AllInclusiveIcon />}
                            disabled={disableSelectAll}
                            onClick={selectAll}
                        />
                    </Tooltip>
                    <Tooltip title="Select None" arrow>
                        <Button
                            variant="contained"
                            color="secondary"
                            className={classes.button}
                            startIcon={<ClearIcon />}
                            disabled={disableSelectNone}
                            onClick={selectNone}
                        />
                    </Tooltip>
                </div>
            );
        } else {
            return(<div></div>); 
        }
    };

    const renderListViewElement = (viewElement, index) => {
        return (<ListViewElement 
                    listType={listType} 
                    key={viewElement.id} 
                    index={index} 
                    id={viewElement.id} 
                    text={viewElement.text} 
                    enabled={viewElement.enabled}
                    setListElementChecked={setListElementChecked}
                    moveListViewElement={moveListViewElement} />);
    };
        
    return (
        <List>
            <ListItem>
                <ListItemText primary={listName} />
            </ListItem>
            <ListItem>
                <ButtonGroup size="medium" aria-label="controls">
                    <Tooltip title={"Submit " + listName + " Changes"} arrow>
                        <Button
                            variant="contained"
                            color="primary"
                            className={classes.button}
                            startIcon={<PlayForWorkIcon />}
                            onClick={handleListViewSubmit}
                        />
                    </Tooltip>
                    {renderGlobalSelect()}
                </ButtonGroup>
            </ListItem>
            <ListItem>
                <FormGroup>
                    {listItemsState.map((element, i) => renderListViewElement(element, i))}
                </FormGroup>
            </ListItem>
        </List>
    );
};

const ListViewElement = ({ listType, id, text, index, enabled, setListElementChecked, moveListViewElement }) => {
    const ref        = useRef(null);

    const handleChange = (event) => {
        setListElementChecked(id, event.target.checked);
    };

    const [{ handlerId }, drop] = useDrop({
        accept: ItemTypes.VIEW_ELEMENT,
        collect(monitor) {
            return {
                handlerId: monitor.getHandlerId(),
            };
        },
        hover(item, monitor) {
            if (!ref.current) {
                return;
            }
            const dragIndex = item.index;
            const hoverIndex = index;
            // Don't replace items with themselves
            if (dragIndex === hoverIndex) {
                return;
            }
            // Determine rectangle on screen
            const hoverBoundingRect = ref.current.getBoundingClientRect();
            // Get vertical middle
            const hoverMiddleY = (hoverBoundingRect.bottom - hoverBoundingRect.top) / 2;
            // Determine mouse position
            const clientOffset = monitor.getClientOffset();
            // Get pixels to the top
            const hoverClientY = clientOffset.y - hoverBoundingRect.top;
            // Only perform the move when the mouse has crossed half of the items height
            // When dragging downwards, only move when the cursor is below 50%
            // When dragging upwards, only move when the cursor is above 50%
            // Dragging downwards
            if (dragIndex < hoverIndex && hoverClientY < hoverMiddleY) {
                return;
            }
            // Dragging upwards
            if (dragIndex > hoverIndex && hoverClientY > hoverMiddleY) {
                return;
            }
            // Time to actually perform the action
            moveListViewElement(dragIndex, hoverIndex);
            // Note: we're mutating the monitor item here!
            // Generally it's better to avoid mutations,
            // but it's good here for the sake of performance
            // to avoid expensive index searches.
            item.index = hoverIndex;
        },
    });
    const [{ isDragging }, drag] = useDrag({
        type: ItemTypes.VIEW_ELEMENT,
        item: () => {
            return { id, index };
        },
        collect: (monitor) => ({
            isDragging: monitor.isDragging(),
        }),
    });
    const opacity = isDragging ? 0 : 1;
    drag(drop(ref));

    return (<FormControlLabel
                ref={ref}
                style={{ opacity }}
                data-handler-id={handlerId}
                control={
                    <Checkbox
                        checked={enabled}
                        onChange={handleChange}
                        name={text}
                        color="primary"
                    />
                }
            label={text} />);
};


const InteractionsSwitch = () => {
    const dispatch = useDispatch();
    const displayInteractions = useSelector(selectDisplayInteractions, shallowEqual);

    const handleChange = (event) => {
        dispatch(setDisplayInteractions(event.target.checked));
    };

    return( <ListItem>
                <FormGroup>
                    <FormControlLabel
                        control={
                                    <Switch
                                        checked={displayInteractions}
                                        onChange={handleChange}
                                        color="primary"
                                    />
                                }
                        label="Display StepwiseQTL Interactions"
                    />
                </FormGroup>
            </ListItem> );
};

const LabelTrackSwitch = () => {
    const dispatch = useDispatch();
    const displayTrackLabels = useSelector(selectDisplayTrackLabels, shallowEqual);

    const handleChange = (event) => {
        dispatch(setDisplayTrackLabels(event.target.checked));
    };

    return( <ListItem>
                <FormGroup>
                    <FormControlLabel
                        control={
                                    <Switch
                                        checked={displayTrackLabels}
                                        onChange={handleChange}
                                        color="primary"
                                    />
                                }
                        label="Display Track Labels"
                    />
                </FormGroup>
            </ListItem> );
};

const ViewController = () => {
    return(
        <List>
            <ConfigurationResetor />
            <Divider />
            <QTLModelCountSelector />
            <Divider />
            <InteractionsSwitch />
            <Divider />
            <LabelTrackSwitch />
            <Divider />
            <QTLConsensusSwitch />
            <Divider />
            <QTLMethodSelector />
            <Divider />
            <ListViewController listName="Linkage Groups" listType="linkageGroups" />
            <Divider />
            <ListViewController listName="Models" listType="models" />
            <Divider />
            <ListViewController listName="Traits"listType="traits" />
        </List>
    );
}; 

export default ViewController;
