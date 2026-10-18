import { useEffect, useId, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../api";
import type { Catalog } from "../api";
import { ErrorMessage, formatCount } from "./Status";

export function TargetCatalog({
  selectedId,
  onSelect,
  catalog,
  compact = false,
}: {
  selectedId: string | null;
  onSelect: (id: string) => void;
  catalog?: Catalog;
  compact?: boolean;
}) {
  const inputId = useId();
  const [input, setInput] = useState("");
  const [filter, setFilter] = useState({
    query: "",
    changed: false,
    offset: 0,
  });
  useEffect(() => {
    const timer = setTimeout(
      () =>
        setFilter((value) =>
          value.query === input.trim()
            ? value
            : { ...value, query: input.trim(), offset: 0 },
        ),
      250,
    );
    return () => clearTimeout(timer);
  }, [input]);
  const { data, error, isFetching } = useQuery({
    queryKey: ["targets", filter],
    queryFn: ({ signal }) =>
      api.targets(filter.query, filter.changed, filter.offset, signal),
  });
  return (
    <section className="catalog" aria-label="Protein catalog">
      {!compact && (
        <div className="catalog-heading">
          <h2>Proteins</h2>
          <span title="Biological records; all are available">
            {catalog ? formatCount(catalog.target_records) : "…"}
          </span>
        </div>
      )}
      <form
        className="catalog-search"
        onSubmit={(event) => {
          event.preventDefault();
          setFilter({ ...filter, query: input.trim(), offset: 0 });
        }}
      >
        <label className="sr-only" htmlFor={inputId}>
          Search gene or transcript
        </label>
        <input
          id={inputId}
          type="search"
          value={input}
          onChange={(e) => setInput(e.target.value)}
          placeholder="Search gene or transcript…"
          autoComplete="off"
        />
        <label
          className="check"
          title="Protein sequence differs from its reference sequence"
        >
          <input
            type="checkbox"
            checked={filter.changed}
            onChange={(e) =>
              setFilter({
                query: input.trim(),
                changed: e.target.checked,
                offset: 0,
              })
            }
          />
          Changed sequences
        </label>
      </form>
      <ErrorMessage error={error} />
      <div className="catalog-count" role="status">
        {isFetching
          ? "Searching…"
          : data
            ? `${formatCount(data.total)} matching records`
            : "Loading…"}
      </div>
      <div
        className="target-list"
        aria-label="Matching proteins"
        aria-busy={isFetching}
      >
        {data?.records.map((row) => (
          <button
            key={row.target_id}
            className="target"
            aria-pressed={selectedId === row.target_id}
            onClick={() => onSelect(row.target_id)}
          >
            <span className="target-name">
              <strong>{row.gene || "Unnamed gene"}</strong>
              <span>{formatCount(row.length)} aa</span>
            </span>
            <span className="target-transcript">{row.transcript}</span>
            <span className="target-origin">
              {row.donor} · {row.target_id.slice(-8)}
              {row.changed ? (
                <span className="changed-dot" title="Changed sequence" />
              ) : null}
            </span>
          </button>
        ))}
        {data?.total === 0 && (
          <p className="muted">
            No matching proteins. Try a gene or transcript prefix.
          </p>
        )}
      </div>
      <div className="pagination">
        <button
          className="quiet"
          disabled={isFetching || !filter.offset}
          onClick={() =>
            setFilter({ ...filter, offset: Math.max(0, filter.offset - 25) })
          }
        >
          Previous
        </button>
        <span>
          {data?.total
            ? `${formatCount(data.offset + 1)}–${formatCount(Math.min(data.offset + data.limit, data.total))}`
            : "0"}
        </span>
        <button
          className="quiet"
          disabled={isFetching || !data || filter.offset + 25 >= data.total}
          onClick={() => setFilter({ ...filter, offset: filter.offset + 25 })}
        >
          Next
        </button>
      </div>
      {catalog && (
        <details className="catalog-details">
          <summary>About this catalog</summary>
          <p>
            {formatCount(catalog.target_records)} biological records ·{" "}
            {formatCount(catalog.unique_sequences)} unique sequences ·{" "}
            {formatCount(catalog.nonreference_records)} changed records.
          </p>
          <p>All records and their provenance are retained.</p>
        </details>
      )}
    </section>
  );
}
