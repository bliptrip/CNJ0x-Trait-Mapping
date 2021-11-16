import React, { useContext, useEffect, useState } from 'react';

import { useSelector } from 'react-redux';

import { makeStyles } from '@material-ui/core/styles';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import {    XGrid, 
            GridToolbar,
            GridApi, 
            useGridApiRef, 
            GridFilterModelParams,
            GridFilterInputValue, 
            GridFilterItem, 
            GridFilterOperator } from '@material-ui/x-grid/dist/index-cjs';

import { selectBlupTableGridFilters } from './blupTableGridSlice';

import blupContext from '../../app/blupContext';

import styles from './BlupTableGrid.module.css';
import DataFrame from 'dataframe-js';

const uuidv4 = require('uuid/v4');

const useStyles = makeStyles({
      table: {
              minWidth: 650,
            },
});

const getBlupTableGridOperators: (t: 'string') => GridFilterOperator[] = (t) => [
  {
    label: '=',
    value: '===',
    getApplyFilterFn: (filterItem: GridFilterItem) => {
      if (!filterItem.columnField || !filterItem.value || !filterItem.operatorValue) {
        return null;
      }

      return ({ value }): boolean => {
        return value === filterItem.value;
      };
    },
    InputComponent: GridFilterInputValue,
    InputComponentProps: { type: t },
  },
  {
    label: '!=',
    value: '!==',
    getApplyFilterFn: (filterItem: GridFilterItem) => {
      if (!filterItem.columnField || !filterItem.value || !filterItem.operatorValue) {
        return null;
      }

      return ({ value }): boolean => {
        return value !== filterItem.value;
      };
    },
    InputComponent: GridFilterInputValue,
    InputComponentProps: { type: t },
  },
  {
    label: '>',
    value: '>',
    getApplyFilterFn: (filterItem: GridFilterItem) => {
      if (!filterItem.columnField || !filterItem.value || !filterItem.operatorValue) {
        return null;
      }

      return ({ value }): boolean => {
        return value > filterItem.value;
      };
    },
    InputComponent: GridFilterInputValue,
    InputComponentProps: { type: t },
  },
  {
    label: '>=',
    value: '>=',
    getApplyFilterFn: (filterItem: GridFilterItem) => {
      if (!filterItem.columnField || !filterItem.value || !filterItem.operatorValue) {
        return null;
      }

      return ({ value }): boolean => {
        return value >= filterItem.value;
      };
    },
    InputComponent: GridFilterInputValue,
    InputComponentProps: { type: t },
  },
  {
    label: '<',
    value: '<',
    getApplyFilterFn: (filterItem: GridFilterItem) => {
      if (!filterItem.columnField || !filterItem.value || !filterItem.operatorValue) {
        return null;
      }

      return ({ value }): boolean => {
        return value < filterItem.value;
      };
    },
    InputComponent: GridFilterInputValue,
    InputComponentProps: { type: t },
  },
  {
    label: '<=',
    value: '<=',
    getApplyFilterFn: (filterItem: GridFilterItem) => {
      if (!filterItem.columnField || !filterItem.value || !filterItem.operatorValue) {
        return null;
      }

      return ({ value }): boolean => {
        return value <= filterItem.value;
      };
    },
    InputComponent: GridFilterInputValue,
    InputComponentProps: { type: t },
  },
];

const BlupStatsTable = ({stats}) => {
    const classes = useStyles();

    return (
        <TableContainer component={Paper}>
            <Table className={classes.table} aria-label="statistics table">
                <TableHead>
                    <TableRow>
                        <TableCell></TableCell>
                        <TableCell align="center">Min</TableCell>
                        <TableCell align="center">Max</TableCell>
                        <TableCell align="center">Range</TableCell>
                        <TableCell align="center">Mean</TableCell>
                        <TableCell align="center">Population SD</TableCell>
                    </TableRow>
                </TableHead>
                <TableBody>
                    {stats.map((s) => (
                        <TableRow key={s.id}>
                            <TableCell component="th" scope="row">
                                {s.id}
                            </TableCell>
                            <TableCell align="right">{s.min}</TableCell>
                            <TableCell align="right">{s.max}</TableCell>
                            <TableCell align="right">{s.range}</TableCell>
                            <TableCell align="right">{s.mean}</TableCell>
                            <TableCell align="right">{s.sdpop}</TableCell>
                        </TableRow>
                    ))}
                </TableBody>
            </Table>
        </TableContainer>
    );
};

export const headCells = [
    { field: 'genotype', type: 'string', headerName: 'Genotype', width: 150, filterOperators: getBlupTableGridOperators('string') },
    { field: 'trait', type: 'string', headerName: 'Trait', width: 150, filterOperators: getBlupTableGridOperators('string') },
    { field: 'model', type: 'string', headerName: 'Model', width: 150, filterOperators: getBlupTableGridOperators('string') },
    { field: 'h2', type: 'number', headerName: 'Genomic Heritability', flex: 1, filterOperators: getBlupTableGridOperators('number') },
    { field: 'vg', type: 'number', headerName: 'Genetic Variance', flex: 1, filterOperators: getBlupTableGridOperators('number') },
    { field: 've', type: 'number', headerName: 'Error Variance', flex: 1, filterOperators: getBlupTableGridOperators('number') },
    { field: 'raw', type: 'number', headerName: 'Value', flex: 1, filterOperators: getBlupTableGridOperators('number') },
    { field: 'blup', type: 'number', headerName: 'BLUP', flex: 1, filterOperators: getBlupTableGridOperators('number') }
];

export default function BlupTableGrid() {
    const apiRef           = useGridApiRef();
    const blupTableData    = useContext(blupContext);
    const filters          = useSelector(selectBlupTableGridFilters);
    const [currFilters, setCurrFilters] = useState([]);
    const [stats,setStats] = useState([]);
    var blupDataDf = undefined;

    useEffect( () => {
        if( apiRef.current.subscribeEvent && blupTableData && blupTableData.length > 0 ) {
            if( blupDataDf === undefined ) {
                blupDataDf = new DataFrame(blupTableData);
            }
            return apiRef.current.subscribeEvent('filterModelChange', (params : GridFilterModelParams) => {
                var filteredData = blupDataDf;
                const filterModel = params.filterModel;
                //Assume all 'And' link operators for now
                filterModel.items.filter( f => f.value ).forEach( f => {
                    filteredData = filteredData.filter( r => eval("r.get(f.columnField) " + f.operatorValue + " f.value")); 
                });
                if( filteredData.toArray().length > 0 ) {
                    const blups = filteredData.stat.stats('blup');
                    const raws = filteredData.stat.stats('raw');
                    setStats( [ {id: 'blup', ...blups, range: blups.max-blups.min}, 
                                {id: 'value', ...raws, range: raws.max-raws.min} ] );
                }
            });
        }
    }, [apiRef, blupTableData] );

    useEffect( () => {
        if( blupTableData && blupTableData.length > 0 && apiRef.current && apiRef.current.upsertFilter ) {
            //Add new filters
            const filtersWithId = filters.map( f => ({ ...f, id: uuidv4() }) );
            filtersWithId.forEach( f => {
                apiRef.current.upsertFilter(f);
            });
            //Delete all previous filters
            currFilters.forEach( f => {
                apiRef.current.deleteFilter(f);
            });
            setCurrFilters(filtersWithId);
        }
    }, [filters, blupTableData] );

    return (
        <div>
            <div>
                <BlupStatsTable stats={stats} />
            </div>
            <div style={{ display: 'flex', height: 400 }}>
                <div style={{ flexGrow: 1 }}>
                    <XGrid
                        apiRef={ apiRef }
                        rows={ blupTableData }
                        columns={ headCells }
                        components={{
                                    Toolbar: GridToolbar,
                                        }}
                    />
                </div>
            </div>
        </div>
    );
}
