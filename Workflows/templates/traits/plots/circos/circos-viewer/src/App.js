/* React and React-Redux libraries */
import React, {useContext, useState, useEffect} from 'react';

/* General support libraries */
import clsx from 'clsx';


/* Material-UI */
import { makeStyles, useTheme } from '@material-ui/core/styles';
import Box from '@material-ui/core/Box';
import CssBaseline from '@material-ui/core/CssBaseline';
import Divider from '@material-ui/core/Divider';
import Drawer from '@material-ui/core/Drawer';
import Grid from '@material-ui/core/Grid';
import Toolbar from '@material-ui/core/Toolbar';
import Tooltip from '@material-ui/core/Tooltip';
import Typography from '@material-ui/core/Typography';

import IconButton from '@material-ui/core/IconButton';
import MenuIcon from '@material-ui/icons/Menu';
import ChevronLeftIcon from '@material-ui/icons/ChevronLeft';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';

/* Application-Specific Features */
import ViewController from './features/viewController/ViewController';
import View from './features/view/View';
import CorrPlot from './features/corrPlot/CorrPlot';
import EffectPlot from './features/effectPlot/EffectPlot';
import LodProfilePlot from './features/lodProfilePlot/LodProfilePlot';
import BlupTable from './features/blupTable/BlupTable';
import BlupTableGrid from './features/blupTableGrid/BlupTableGrid';

import * as d3 from 'd3';
import {DataFrame} from 'dataframe-js';

import blupContext from './app/blupContext';

import './App.css';

const drawerWidth = 300;

const useStyles = makeStyles((theme) => ({
  root: {
    display: 'flex',
  },
  appBar: {
    borderRadius: 10,
    backgroundColor: "#EEEFFFD0",
    transition: theme.transitions.create(['margin', 'width'], {
      easing: theme.transitions.easing.sharp,
      duration: theme.transitions.duration.leavingScreen,
    }),
  },
  menuButton: {
    position: "relative",
    backgroundColor: "#FFFFFFC0",
    left: 13,
    top: 13,
    marginRight: theme.spacing(2),
    zIndex: 2000
  },
  hide: {
    visibility: 'hidden',
  },
  drawer: {
    width: drawerWidth,
    flexShrink: 0,
  },
  drawerPaper: {
    width: drawerWidth,
  },
  drawerHeader: {
    display: 'flex',
    alignItems: 'center',
    padding: theme.spacing(0, 1),
    // necessary for content to be below app bar
    ...theme.mixins.toolbar,
    justifyContent: 'flex-end',
  },
  content: {
    width: "90%",
    padding: theme.spacing(0),
    transition: theme.transitions.create('margin', {
      easing: theme.transitions.easing.sharp,
      duration: theme.transitions.duration.leavingScreen,
    }),
    marginLeft: -drawerWidth,
  },
  contentShift: {
    width: `calc(100% - ${drawerWidth}px)`,
    transition: theme.transitions.create('margin', {
      easing: theme.transitions.easing.easeOut,
      duration: theme.transitions.duration.enteringScreen,
    }),
    marginLeft: 0,
  },
}));

function PersistentDrawerLeft() {
  const classes = useStyles();
  const theme = useTheme();
  const [open, setOpen] = React.useState(false);

  const handleDrawerOpen = () => {
    setOpen(true);
  };

  const handleDrawerClose = () => {
    setOpen(false);
  };

  return (
    <div className={classes.root}>
      <CssBaseline />
      <Box
        className={clsx(classes.appBar, {
          [classes.hide]: open,
        })}
      >
        <Tooltip title="Configure Graph" arrow>
            <IconButton
                color="inherit"
                aria-label="open drawer"
                onClick={handleDrawerOpen}
                edge="start"
                className={clsx(classes.menuButton, open && classes.hide)}
            >
                <MenuIcon />
            </IconButton>
        </Tooltip>
      </Box>
      <Drawer
        className={classes.drawer}
        anchor="left"
        variant="persistent"
        open={open}
        onClose={handleDrawerClose}
        classes={{
          paper: classes.drawerPaper,
        }}
      >
        <div className={classes.drawerHeader}>
          <IconButton onClick={handleDrawerClose}>
            {theme.direction === 'ltr' ? <ChevronLeftIcon /> : <ChevronRightIcon />}
          </IconButton>
        </div>
        <Divider />
        <ViewController />
      </Drawer>
      <main
        className={clsx(classes.content)}
        onClick={handleDrawerClose}
      >
        <Grid 
            container 
            spacing={1}
            direction="row"
            justify="space-between"
            alignItems="flex-start">
            <Grid item xs={9}>
                <View />
            </Grid>
            <Grid 
                container 
                item
                xs={3}
                spacing={1}
                direction="column"
                justify="flex-start"
                alignItems="center">
                <Grid item xs={12}>
                    <EffectPlot />
                </Grid>
                <Grid item xs={12}>
                    <LodProfilePlot />
                </Grid>
                <Grid item xs={12}>
                    <CorrPlot />
                </Grid>
            </Grid>
            <Grid item xs={12}>
                <BlupTableGrid />
            </Grid>
        </Grid>
      </main>
    </div>
  );
}

function App() {
    const [init, setInit]     = useState(false);
    const [blups, setBlups]   = useState([]);

    useEffect( () => {
        d3.json('configs/blups_collated.wide.json')
        .then( d => {
            setBlups(d);
        });
    }, [init]);

    if( init === false ) {
        setInit(true);
    }

    return (
        <div className="App">
            <blupContext.Provider value={blups}>
                <PersistentDrawerLeft />
            </blupContext.Provider>
        </div>
    );
}

export default App;
