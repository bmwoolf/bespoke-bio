import { useEffect, useId, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../api";
import { ErrorMessage, formatCount } from "./Status";

export function DesignCatalog({
  selectedId,
  onSelect,
}: {
  selectedId: string | null;
  onSelect: (id: string) => void;
}) {
  const inputId = useId();
  const [input, setInput] = useState("");
  const [filter, setFilter] = useState({ query: "", offset: 0 });
  useEffect(() => {
    const timer = setTimeout(
      () =>
        setFilter((value) =>
          value.query === input.trim()
            ? value
            : { query: input.trim(), offset: 0 },
        ),
      250,
    );
    return () => clearTimeout(timer);
  }, [input]);
  const { data, error, isFetching } = useQuery({
    queryKey: ["designs", filter.query, filter.offset],
    queryFn: ({ signal }) => api.designs(filter.query, filter.offset, signal),
  });
  return (
    <section className="catalog design-catalog" aria-label="Designed targets">
      <div className="catalog-search">
        <label className="sr-only" htmlFor={inputId}>
          Search designed targets
        </label>
        <input
          id={inputId}
          type="search"
          placeholder="Search target or variant…"
          value={input}
          onChange={(event) => setInput(event.target.value)}
          autoComplete="off"
        />
      </div>
      <ErrorMessage error={error} />
      <div className="catalog-count" role="status">
        {isFetching
          ? "Loading…"
          : `${formatCount(data?.total ?? 0)} saved jobs`}
      </div>
      <div
        className="target-list"
        aria-label="Matching designs"
        aria-busy={isFetching}
      >
        {data?.records.map((row) => (
          <button
            key={row.job_id}
            className="target"
            aria-pressed={selectedId === row.job_id}
            onClick={() => onSelect(row.job_id)}
          >
            <span className="target-name">
              <strong>{row.name}</strong>
              <span>{row.candidate_count} keys</span>
            </span>
            <span className="target-transcript">{row.description}</span>
            <span className="target-availability">
              {row.length ? `${row.length} aa · ` : ""}
              {row.passed_count} pass native filters
            </span>
            <span className="target-origin">
              {row.source} · {row.job_id.slice(0, 6)}
            </span>
          </button>
        ))}
        {data?.total === 0 && (
          <p className="muted">
            {filter.query
              ? "No matching designs."
              : "Completed batch results will appear here once imported."}
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
            ? `${data.offset + 1}–${Math.min(data.offset + data.limit, data.total)}`
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
      <p className="catalog-footnote">
        Saved predictions. Filter passes do not establish binding.
      </p>
    </section>
  );
}
