import React, { useEffect, useRef, useState, useCallback } from 'react';
import { useSelector, shallowEqual, useDispatch } from 'react-redux';
import { useDrag, useDrop } from 'react-dnd';
import {selectList,
        selectListItem,
        selectDisplayTrackLabels,
        selectQTLConsensus,
        selectQTLMethod,
        setQTLConsensus,
        setQTLMethod,
        setList,
        setListItemChecked,
        setDisplayTrackLabels
} from './viewControllerSlice';
//import styles from './ViewController.module.css';
import update from 'immutability-helper';

/* Material-UI Core and Components */
import Button from '@material-ui/core/Button';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import Checkbox from '@material-ui/core/Checkbox';
import Divider from '@material-ui/core/Divider';
import FormGroup from '@material-ui/core/FormGroup';
import Radio from '@material-ui/core/Radio';
import RadioGroup from '@material-ui/core/RadioGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import FormControl from '@material-ui/core/FormControl';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import Switch from '@material-ui/core/Switch';
import Tooltip from '@material-ui/core/Tooltip';

import { makeStyles } from '@material-ui/core/styles';

/* Material-UI Icons */
import AllInclusiveIcon from '@material-ui/icons/AllInclusive';
import ClearIcon from '@material-ui/icons/Clear';
import PlayForWorkIcon from '@material-ui/icons/PlayForWork';

const ItemTypes = {
    VIEW_ELEMENT: 'ViewElement',
};

const useStyles = makeStyles((theme) => ({
    button: {
        margin: theme.spacing(0.5),
    },
}));

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
                <React.Fragment>
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
                </React.Fragment>
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
